from netvis import *
import pandas as pd
import numpy as np
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




cmp_df_selected_with_index = pd.read_csv('data/cmp_df_selected_with_index.csv')
org_df_selected_with_index = pd.read_csv('data/org_df_selected_with_index.csv')
org_cmp_selected_metapath = pd.read_csv('data/org_cmp_manually_selected_metapath.csv')
dwpc_mat_2d = np.load('data/dwpc_mat_2d.npy')

compound_names = list(cmp_df_selected_with_index.compound_name.unique())
compound_names.sort()

MAX_COUNT = org_df_selected_with_index.shape[0]


def main():    
    st.markdown("<h1 style='text-align: center; color: black;'>BCMM Compounds - SPOKE insight</h1>", unsafe_allow_html=True)
    st.sidebar.header("Sample Compounds")
    hyperlink_names = ["Cholic acid", "ursodiol", "estrone 3-sulfate", "doxorubicin", "lithocholic acid"]
    compound_selected_sample = st.sidebar.radio("", [None] + hyperlink_names)
    st.sidebar.header("Full Search - Select Compound")
    compound_selected_search = get_search_term(compound_names)
    bacteria_count = sidebar_options()
    if compound_selected_search == DEFAULT_SELECTION and compound_selected_sample != None:
        write_bacteria_table(compound_selected_sample, bacteria_count)
        st.markdown("<h4 style='text-align: left; color: black;'>Explore SPOKE network for {}</h4>".format(compound_selected_sample), unsafe_allow_html=True)
        print_vpn_warning()
        organism_id = st.text_input("Enter NCBI ID of the Organism")
        compound_id = cmp_df_selected_with_index[cmp_df_selected_with_index.compound_name==compound_selected_sample].spoke_identifier.values[0]
        metapath_based_network_vis(organism_id, compound_id, org_cmp_selected_metapath)
        
        
    if compound_selected_search != DEFAULT_SELECTION and compound_selected_sample != None:
        st.markdown("<h5 style='text-align: center; color: black;'>Multiple Compound selection was made.</h5>", unsafe_allow_html=True)
        st.markdown("<h5 style='text-align: center; color: black;'>Reset either Sample or Search Compound to None</h5>", unsafe_allow_html=True)
    if compound_selected_search != DEFAULT_SELECTION and compound_selected_sample == None:
        write_bacteria_table(compound_selected_search, bacteria_count)
        st.markdown("<h4 style='text-align: left; color: black;'>Explore SPOKE network for {}</h4>".format(compound_selected_search), unsafe_allow_html=True)
        print_vpn_warning()
        organism_id = st.text_input("Enter NCBI ID of the Organism")
        compound_id = cmp_df_selected_with_index[cmp_df_selected_with_index.compound_name==compound_selected_search].spoke_identifier.values[0]
        metapath_based_network_vis(organism_id, compound_id, org_cmp_selected_metapath)
        
        
    
    
def sidebar_options():
    bacteria_count = st.sidebar.slider('Bacteria count', MIN_COUNT, MAX_COUNT, DEFAULT_COUNT)
    return bacteria_count
    

def write_bacteria_table(compound_selected, bacteria_count):    
    st.markdown("<h4 style='text-align: left; color: black;'>Top Bacterial nodes associated with {}</h4>".format(compound_selected), unsafe_allow_html=True)
    st.write(get_bacteria_table(compound_selected, bacteria_count))
    
    
def get_bacteria_table(cmp_name, bacteria_count):
    column_ind = cmp_df_selected_with_index[cmp_df_selected_with_index.compound_name==cmp_name].cmp_index.values[0]
    cmp_df_selected_with_index_copy = org_df_selected_with_index.copy()
    cmp_df_selected_with_index_copy.loc[:,'dwpc'] = dwpc_mat_2d[:,column_ind]
    cmp_df_selected_with_index_copy = cmp_df_selected_with_index_copy[['spoke_id', 'spoke_name', 'dwpc']]
    cmp_df_selected_with_index_copy.spoke_id = cmp_df_selected_with_index_copy.spoke_id.astype(str)   
    cmp_df_selected_with_index_copy = cmp_df_selected_with_index_copy.rename(columns={'spoke_id':'NCBI ID', 'spoke_name':'name'})
    cmp_df_selected_with_index_copy = cmp_df_selected_with_index_copy.sort_values(by='dwpc', ascending=False).head(bacteria_count)
    return cmp_df_selected_with_index_copy.reset_index().drop("index", axis=1)
    
    
def get_search_term(compound_names):
    return st.sidebar.selectbox(" ", [DEFAULT_SELECTION]+compound_names, index=0)


if __name__ == "__main__":
    main()
