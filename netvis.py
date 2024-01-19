import streamlit as st
import streamlit.components.v1 as components
from neo4j import GraphDatabase, basic_auth
import networkx as nx
from pyvis.network import Network
import os
import tempfile
from dotenv import load_dotenv


node_color_map = {
        "Anatomy": "#b3de69",
        "BiologicalProcess": "#fdb462",
        "Blend": "#fcd9fc",
        "CellType": "#6b853f",
        "CellularComponent": "#ffffb3",
        "ClinicalLab": "#e2c3b8",
        "Complex": "#4575b4",
        "Compound": "#bc80bd",
        "Disease": "#fb8072",
        "EC": "#f79767",
        "Food": "#c7b78f",
        "Gene": "#80b1d3",
        "Location": "#e17910",
        "MolecularFunction": "#ffed6f",
        "MiRNA": "#f3ecc5",
        "Organism": "#d9c8ae",
        "Pathway": "#57c7e3",
        "PharmacologicClass": "#bebada",
        "Protein": "#8dd3c7",
        "ProteinDomain": "#00ffff",
        "ProteinFamily": "#006666",
        "PwGroup": "#5D3FD3",
        "Reaction": "#f16667",
        "SideEffect": "#ccebc5",
        "Symptom": "#fccde5"
    }


load_dotenv(os.path.join(os.path.expanduser("~"), ".neo4j_config.env"))
URL = os.environ.get("SPOKE_URI")
SPOKE_USER = os.environ.get("SPOKE_USER")
SPOKE_PASSWORD = os.environ.get("SPOKE_PSW")

def connect_to_neo4j():
    return GraphDatabase.driver(URL, auth=basic_auth(SPOKE_USER, SPOKE_PASSWORD))  

def fetch_shortest_path(driver, source, target):
    with driver.session() as session:
        result = session.run(
            """            
            MATCH (o:Organism {identifier: $source}), (c:Compound {identifier: $target})
            MATCH path = allShortestPaths((c)-[*]-(o))
            RETURN path
            """,
            source=source,
            target=target,
        )
        paths = []
        for row in result:
            paths.append(row["path"])
        return paths

def fetch_path(driver, source, target):
    with driver.session() as session:
        result = session.run(
            """            
            MATCH (o:Organism {identifier: $source}), (c:Compound {identifier: $target})
            MATCH path1 = (o)-[:ENCODES_OeP]->(p:Protein)-[:HAS_PhEC]->(e:EC)-[:ISA_ECiEC]->(e2:EC)-[*1..2]->(r:Reaction)
            MATCH path2 = allShortestPaths((r)-[*]-(c))            
            RETURN path1, path2 LIMIT 10
            """,
            source=source,
            target=target,
        )
        paths = []
        for row in result:
            paths.append(row["path1"])
            paths.append(row["path2"])
        result = session.run(
            """            
            MATCH (o:Organism {identifier: $source}), (c:Compound {name: $target})
            MATCH path1 = (o)-[:ENCODES_OeP]->(p:Protein)-[:HAS_PhEC]->(e:EC)-[*1..2]->(r:Reaction)
            MATCH path2 = allShortestPaths((r)-[*]-(c))            
            RETURN path1, path2 LIMIT 10
            """,
            source=source,
            target=target,
        )        
        for row in result:
            paths.append(row["path1"])
            paths.append(row["path2"])
        if len(paths) == 0:
            paths = fetch_shortest_path(driver, source, target)
        return paths
    
    
    
def fetch_path_from_metapath(driver, source_node, target_node, metapath_data):
    nhop = metapath_data['nhops']
    column_names = [f"col_{i}" for i in range(1, nhop+1)]
    place_holder_values = (source_node,) + tuple(metapath_data[col] for col in column_names) + (target_node,)    
    cypher = 'MATCH path=(o:Organism{identifier:%s})-[:%s]->(n1)'
    for i in range(nhop-2):
        cypher += '-[:%s]-(n{})'.format(i+2)
    cypher += '-[:%s]-(c:Compound{identifier:"%s"})'
    cypher += ' UNWIND nodes(path) AS n'
    cypher += ' UNWIND relationships(path) AS r'
    cypher += ' RETURN n, r'
    cypher = cypher%(place_holder_values)
    with driver.session() as session:
        with session.begin_transaction() as tx:
            result = tx.run(cypher)
            path_list = []
            for row in result:
                path_list.append((row['n'], row['r']))
    return path_list


    
    
def create_legend(legend_color_map):
    legend_html = """
    <div style="border: 2px solid #ccc; padding: 10px; background-color: #f9f9f9;">
        <h5>Color Legend</h5>
        <div style="display: flex; flex-direction: column; align-items: flex-start;">
    """

    for nodetype, color in legend_color_map.items():
        legend_html += f"""
            <div style="display: flex; align-items: center;">
                <div style="width: 20px; height: 20px; background-color: {color};"></div>
                <span style="margin-left: 5px;">{nodetype}</span>
            </div>
        """

    legend_html += """
        </div>
    </div>
    """

    return legend_html

def create_nx_graph(paths):
    graph = nx.DiGraph()
    for path in paths:
        for record in path:            
            sub_obj = []
            for node in record.nodes:
                node_label = list(node.labels)[0]
                try:
                    graph.add_node(node["name"], nodetype=node_label)
                    sub_obj.append(node["name"])
                except:
                    graph.add_node(node["identifier"], nodetype=node_label)
                    sub_obj.append(node["identifier"])
            graph.add_edge(sub_obj[0], sub_obj[1], edgetype=record.type)
    return graph

def create_nx_graph_v2(paths):
    graph = nx.DiGraph()
    for path in paths:
        path = list(map(lambda x:x[-1], path))
        path = list(set(path))
        for record in path:            
            sub_obj = []
            for node in record.nodes:
                node_label = list(node.labels)[0]
                try:
                    if node_label != 'Protein':
                        graph.add_node(node["name"], nodetype=node_label)
                        sub_obj.append(node["name"])
                    else:
                        graph.add_node(node["description"], nodetype=node_label)
                        sub_obj.append(node["description"])
                except:
                    graph.add_node(node["identifier"], nodetype=node_label)
                    sub_obj.append(node["identifier"])
            graph.add_edge(sub_obj[0], sub_obj[1], edgetype=record.type)
    return graph



def create_pyvis_graph(graph):
    net = Network(width="100%", height="800px", notebook=True)
    net.from_nx(graph)
    net_color_mapped = []
    legend_color_map = {}
    for item in net.nodes:
        item["color"] = node_color_map[item["nodetype"]]
        net_color_mapped.append(item)
        legend_color_map[item["nodetype"]] = item["color"]
    net.nodes = net_color_mapped
    for edge in net.edges:
        edge["arrows"] = "to"
        source_node = edge["from"]
        target_node = edge["to"]
        edge_type = graph.edges[source_node, target_node]["edgetype"]
        edge["title"] = f"Edge Type: {edge_type}"
    net.repulsion(
                    node_distance=250,
                    central_gravity=0.33,
                    spring_length=110,
                    spring_strength=0.1,
                    damping=100
                   )
    return net, legend_color_map

def show_network(net, legend_color_map):                               
    with tempfile.NamedTemporaryFile(delete=False, suffix=".html") as temp_html_file:
        temp_html_filename = temp_html_file.name
        net.save_graph(temp_html_filename)
    legend_html = create_legend(legend_color_map)
    with open(temp_html_filename, 'r', encoding='utf-8') as temp_html_file:
        components.html(legend_html + temp_html_file.read(), height=1000, width=1000)
    os.unlink(temp_html_filename)

    
def network_vis(organism_id, compound_id):
    if organism_id and compound_id: 
        with st.spinner("Connecting to SPOKE ..."):
            driver = connect_to_neo4j()
            paths = fetch_path(driver, int(organism_id), compound_id)
            driver.close()
            graph = create_nx_graph(paths)
            net, legend_color_map = create_pyvis_graph(graph)
            show_network(net, legend_color_map)

def metapath_based_network_vis(organism_id, compound_id, metapath_data):
    if organism_id and compound_id: 
        with st.spinner("Connecting to SPOKE ..."):
            driver = connect_to_neo4j()
            paths = []
            for index, row in metapath_data.iterrows():
                paths.append(fetch_path_from_metapath(driver, int(organism_id), compound_id, row))
            driver.close()
            graph = create_nx_graph_v2(paths)
            net, legend_color_map = create_pyvis_graph(graph)
            show_network(net, legend_color_map)
            
            
def print_vpn_warning():
    st.markdown(
    '<div style="color: red;">IMPORTANT : Make sure you are connected to UCSF VPN to explore the network</div>',
    unsafe_allow_html=True
)
    