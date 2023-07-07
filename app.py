import streamlit as st
import pandas as pd
import pickle


DEFAULT_SELECTION = "None"
DEFAULT_COUNT = 25
MIN_COUNT = DEFAULT_COUNT


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
    data_selected = data[compound_selected]
    data_selected["ncbi_id"] = data["ncbi_id"]
    data_selected["name"] = data["name"]
    data_selected_df = pd.DataFrame(data_selected)
    data_selected_df = data_selected_df[["ncbi_id", "name", "embedding", "p_value"]]
    if sort_by == "embedding score":
        bacteria_df = data_selected_df.sort_values(by="embedding", ascending=False).head(bacteria_count)
    else:
        bacteria_df = data_selected_df.sort_values(by="p_value", ascending=True).head(bacteria_count)
        
    bacteria_df.ncbi_id = bacteria_df.ncbi_id.astype(str)    
    bacteria_df.rename(columns={"ncbi_id": "NCBI ID", "embedding": "embedding score", "p_value": "proximity pvalue"}, inplace=True)                
    return bacteria_df.reset_index().drop("index", axis=1)
    

def get_search_term(compound_names):
    return st.sidebar.selectbox(" ", [DEFAULT_SELECTION]+compound_names, index=0)
        



if __name__ == "__main__":
    main()
