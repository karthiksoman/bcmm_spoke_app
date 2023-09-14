from netvis import *
import pandas as pd
import pickle
import plotly.express as px
import plotly.graph_objects as go


DEFAULT_SELECTION = "None"
DEFAULT_COUNT = 25
MIN_COUNT = DEFAULT_COUNT
FIG_OPACITY = 1
MARKER_SIZE = 5
FIG_WIDTH = 500
FIG_HEIGHT = 700
FONT_SIZE = 18
LEGEND_SIZE = 12


with open("data/bcmm_compounds_all_bacteria_with_proximity_pvalue.pickle", "rb") as f:
    data = pickle.load(f)
compound_names = list(data.keys())[:-2]
compounds_to_remove = ['diacetamate',
 'azetirelin',
 'ipsalazide',
 'indicine n-oxide',
 'cyclohexylsulfamate',
 'cortisol']
compound_names = list(set(compound_names) - set(compounds_to_remove))
compound_names.sort()
MAX_COUNT = data[compound_names[0]]["embedding"].shape[0]

cmp_map = pd.read_csv("data/bcmm_compounds_combined_refined.csv")


def main():    
    st.markdown("<h1 style='text-align: center; color: black;'>BCMM Compounds - SPOKE insight</h1>", unsafe_allow_html=True)
    st.sidebar.header("Sample Compounds")
    hyperlink_names = ["dolasetron", "fludrocortisone acetate", "phenazopyridine", "trandolapril", "biperiden"]
    compound_selected_sample = st.sidebar.radio("", [None] + hyperlink_names)
    st.sidebar.header("Full Search - Select Compound")
    compound_selected_search = get_search_term(compound_names)    
    bacteria_count = sidebar_options()
    if compound_selected_search == DEFAULT_SELECTION and compound_selected_sample != None:
        write_bacteria_table(compound_selected_sample, bacteria_count)
        st.markdown("<h4 style='text-align: left; color: black;'>Explore SPOKE network for {}</h4>".format(compound_selected_sample), unsafe_allow_html=True)
        print_vpn_warning()
        organism_id = st.text_input("Enter NCBI ID of the Organism")
        compound_id = cmp_map[cmp_map.compound_name==compound_selected_sample].spoke_identifier.values[0]
        network_vis(organism_id, compound_id)

    if compound_selected_search != DEFAULT_SELECTION and compound_selected_sample != None:
        st.markdown("<h5 style='text-align: center; color: black;'>Multiple Compound selection was made.</h5>", unsafe_allow_html=True)
        st.markdown("<h5 style='text-align: center; color: black;'>Reset either Sample or Search Compound to None</h5>", unsafe_allow_html=True)
    if compound_selected_search != DEFAULT_SELECTION and compound_selected_sample == None:
        write_bacteria_table(compound_selected_search, bacteria_count)
        st.markdown("<h4 style='text-align: left; color: black;'>Explore SPOKE network for {}</h4>".format(compound_selected_search), unsafe_allow_html=True)
        print_vpn_warning()
        organism_id = st.text_input("Enter NCBI ID of the Organism")
        compound_id = cmp_map[cmp_map.compound_name==compound_selected_search].spoke_identifier.values[0]
        network_vis(organism_id, compound_id)

def print_vpn_warning():
    st.markdown(
    '<div style="color: red;">IMPORTANT : Make sure you are connected to UCSF VPN to explore the network</div>',
    unsafe_allow_html=True
)


def sidebar_options():
        bacteria_count = st.sidebar.slider('Bacteria count', MIN_COUNT, MAX_COUNT, DEFAULT_COUNT)
        st.sidebar.markdown("<h4 style='text-align: left; color: black;'>What is embedding score?</h4>", unsafe_allow_html=True)
        st.sidebar.markdown("<h5 style='text-align: left; color: black;'>It is the value obtained from applying personalized page rank algorithm on SPOKE graph. It represents the relative saliency of a bacterial node to the selected compound (expressed as percentile score)</h5>", unsafe_allow_html=True)
        return bacteria_count
    
def write_bacteria_table(compound_selected, bacteria_count):    
    st.markdown("<h4 style='text-align: left; color: black;'>Top Bacterial nodes associated with {}</h4>".format(compound_selected), unsafe_allow_html=True)
    st.write(get_bacteria_table(compound_selected, bacteria_count))

    
    
    
                 
def get_bacteria_table(compound_selected, bacteria_count):
    data_selected = data[compound_selected]
    data_selected["ncbi_id"] = data["ncbi_id"]
    data_selected["name"] = data["name"]
    data_selected_df = pd.DataFrame(data_selected)
    data_selected_df = data_selected_df[["ncbi_id", "name", "embedding"]]
    bacteria_df = data_selected_df.sort_values(by="embedding", ascending=False).head(bacteria_count)        
    bacteria_df.ncbi_id = bacteria_df.ncbi_id.astype(str)    
    bacteria_df.rename(columns={"ncbi_id": "NCBI ID", "embedding": "embedding score"}, inplace=True)                
    return bacteria_df.reset_index().drop("index", axis=1)
    

def get_search_term(compound_names):
    return st.sidebar.selectbox(" ", [DEFAULT_SELECTION]+compound_names, index=0)
        



if __name__ == "__main__":
    main()
