import streamlit as st
import pandas as pd
import json


DEFAULT_SELECTION = "None"
MIN_COUNT = 5
MAX_COUNT = 50
DEFAULT_COUNT = 10


with open("data/bcmm_compounds_top_bacteria_with_proximity_score.json", "r") as f:
    json_data = json.load(f)
compound_names = list(json_data.keys())
compounds_to_remove = ['diacetamate',
 'azetirelin',
 'ipsalazide',
 'indicine n-oxide',
 'cyclohexylsulfamate',
 'cortisol']
compound_names = list(set(compound_names) - set(compounds_to_remove))
compound_names.sort()


def main():    
    st.markdown("<h1 style='text-align: center; color: black;'>BCMM Compounds - SPOKE insight</h1>", unsafe_allow_html=True)
    st.sidebar.header("Sample Compounds")
    hyperlink_names = ["dolasetron", "fludrocortisone acetate", "phenazopyridine", "trandolapril", "biperiden"]
    compound_selected_sample = st.sidebar.radio("", [None] + hyperlink_names)
    st.sidebar.header("Full Search - Select Compound")
    compound_selected_search = get_search_term(compound_names)    
    bacteria_count, sort_by = sidebar_options()
    if compound_selected_search == DEFAULT_SELECTION and compound_selected_sample != None:
        write_bacteria_table(compound_selected_sample, bacteria_count, sort_by)
    if compound_selected_search != DEFAULT_SELECTION and compound_selected_sample != None:
        st.markdown("<h5 style='text-align: center; color: black;'>Multiple Compound selection was made.</h5>", unsafe_allow_html=True)
        st.markdown("<h5 style='text-align: center; color: black;'>Reset either Sample or Search Compound to None</h5>", unsafe_allow_html=True)
    if compound_selected_search != DEFAULT_SELECTION and compound_selected_sample == None:
        write_bacteria_table(compound_selected_search, bacteria_count, sort_by)



def sidebar_options():
        bacteria_count = st.sidebar.slider('Bacteria count', MIN_COUNT, MAX_COUNT, DEFAULT_COUNT)
        sort_by = st.sidebar.selectbox("How to sort", ["embedding score", "proximity pvalue"], index=1)
        st.sidebar.markdown("<h5 style='text-align: left; color: black;'>embedding score: It represents the relative saliency of a bacterial node to the selected compound (expressed as percentile score)</h5>", unsafe_allow_html=True)
        st.sidebar.markdown("<h5 style='text-align: left; color: black;'>proximity pvalue: It represents how proximal is the bacterial node to the selected compound in SPOKE graph, compared to a random bacterial node</h5>", unsafe_allow_html=True)
        return bacteria_count, sort_by
    
def write_bacteria_table(compound_selected, bacteria_count, sort_by):    
    st.markdown("<h4 style='text-align: left; color: black;'>Top Bacterial nodes associated with {}</h4>".format(compound_selected), unsafe_allow_html=True)
    st.write(get_bacteria_table(compound_selected, bacteria_count, sort_by))

                        
                 
def get_bacteria_table(compound_selected, bacteria_count, sort_by):
    if sort_by == "embedding score":
        bacteria_df = pd.DataFrame(json_data[compound_selected]).sort_values(by="percentile_score", ascending=False).head(bacteria_count)
    else:
        bacteria_df = pd.DataFrame(json_data[compound_selected]).sort_values(by="proximity_pvalue", ascending=True).head(bacteria_count)
    bacteria_df.ncbi_id = bacteria_df.ncbi_id.astype(str)    
    bacteria_df.rename(columns={"ncbi_id": "NCBI ID", "percentile_score": "embedding score", "proximity_pvalue": "proximity pvalue"}, inplace=True)                
    return bacteria_df.reset_index().drop("index", axis=1)
    

def get_search_term(compound_names):
    return st.sidebar.selectbox(" ", [DEFAULT_SELECTION]+compound_names, index=0)
        



if __name__ == "__main__":
    main()
