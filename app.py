import streamlit as st
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
        plot_bacteria_table(compound_selected_sample)
    if compound_selected_search != DEFAULT_SELECTION and compound_selected_sample != None:
        st.markdown("<h5 style='text-align: center; color: black;'>Multiple Compound selection was made.</h5>", unsafe_allow_html=True)
        st.markdown("<h5 style='text-align: center; color: black;'>Reset either Sample or Search Compound to None</h5>", unsafe_allow_html=True)
    if compound_selected_search != DEFAULT_SELECTION and compound_selected_sample == None:
        write_bacteria_table(compound_selected_search, bacteria_count, sort_by)
        plot_bacteria_table(compound_selected_search)



def sidebar_options():
        bacteria_count = st.sidebar.slider('Bacteria count', MIN_COUNT, MAX_COUNT, DEFAULT_COUNT)
        sort_by = st.sidebar.selectbox("How to sort", ["embedding score", "proximity pvalue"], index=1)
        st.sidebar.markdown("<h5 style='text-align: left; color: black;'>embedding score: It represents the relative saliency of a bacterial node to the selected compound (expressed as percentile score)</h5>", unsafe_allow_html=True)
        st.sidebar.markdown("<h5 style='text-align: left; color: black;'>proximity pvalue: It represents how proximal is the bacterial node to the selected compound in SPOKE graph, compared to a random bacterial node</h5>", unsafe_allow_html=True)
        return bacteria_count, sort_by
    
def write_bacteria_table(compound_selected, bacteria_count, sort_by):    
    st.markdown("<h4 style='text-align: left; color: black;'>Top Bacterial nodes associated with {}</h4>".format(compound_selected), unsafe_allow_html=True)
    st.write(get_bacteria_table(compound_selected, bacteria_count, sort_by))


def plot_bacteria_table(compound_selected):
    data_selected = data[compound_selected]
    data_selected["ncbi_id"] = data["ncbi_id"]
    data_selected["name"] = data["name"]
    data_selected_plot = pd.DataFrame(data_selected)
    data_selected_plot = data_selected_plot[["ncbi_id", "name", "embedding", "p_value"]]
    data_selected_plot["color"] = "gray"
    data_selected_plot.loc[data_selected_plot["p_value"] < 0.05, "color"] = "red"
    data_selected_plot["size"] = MARKER_SIZE
    color_discrete_map = {"gray": "gray", "red": "red"}

    fig = px.scatter(data_selected_plot, 
                     x="p_value", y="embedding", 
                     color="color",
                     color_discrete_map=color_discrete_map,
                     opacity=FIG_OPACITY, 
                     hover_name="name",
                     hover_data=["name"]
                    )

    fig.update_layout(
        margin=dict(l=20, r=20, t=20, b=20),
        width=FIG_WIDTH + 100,
        height=FIG_HEIGHT - 300,
        showlegend=True,
        legend=dict(
            title="Legend",
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1,
            font=dict(
                size=LEGEND_SIZE
            )
        ),
        xaxis=dict(
            showgrid=True,
            showticklabels=True,
            tickfont=dict(size=FONT_SIZE),
            title=dict(
                text="Proximity p-value",
                font=dict(size=FONT_SIZE)
            )
        ),
        yaxis=dict(
            showgrid=True,
            showticklabels=True,
            tickfont=dict(size=FONT_SIZE),
            title=dict(
                text="Embedding score",
                font=dict(size=FONT_SIZE)
            )
        )
    )

    fig.update_traces(marker=dict(size=MARKER_SIZE))
    fig.update_traces(showlegend=False)

    # Add desired legend entries
    fig.add_trace(go.Scatter(
        x=[None],
        y=[None],
        mode="markers",
        marker=dict(size=MARKER_SIZE, color="red"),
        name="Bacteria significantly proximal to {} in SPOKE graph".format(compound_selected)
    ))
    fig.add_trace(go.Scatter(
        x=[None],
        y=[None],
        mode="markers",
        marker=dict(size=MARKER_SIZE, color="gray"),
        name="Bacteria NOT significantly proximal to {} in SPOKE graph".format(compound_selected)
    ))
    st.markdown("<h4 style='text-align: center; color: black;'>Bacterial distribution in embedding and p-value space (associated with {})</h4>".format(compound_selected), unsafe_allow_html=True)
    st.plotly_chart(fig)
    
    
# def plot_bacteria_table(compound_selected):
#     data_selected = data[compound_selected]
#     data_selected["ncbi_id"] = data["ncbi_id"]
#     data_selected["name"] = data["name"]
#     data_selected_plot = pd.DataFrame(data_selected)
#     data_selected_plot = data_selected_plot[["ncbi_id", "name", "embedding", "p_value"]]
#     data_selected_plot["color"] = "gray"
#     data_selected_plot.loc[data_selected_plot["p_value"] < 0.05, "color"] = "red"
#     data_selected_plot["size"] = MARKER_SIZE
#     color_discrete_map = {"gray": "gray", "red": "red"}
#     fig = px.scatter(data_selected_plot, 
#                  x="p_value", y="embedding", 
#                  color="color",
#                  color_discrete_map=color_discrete_map,
#                  opacity=FIG_OPACITY, 
#                  hover_name="name",
#                  hover_data=["name"]
#                 )
#     fig_final = go.Figure(data=fig.data)
#     fig_final.update_layout(
#                             margin=dict(l=20, r=20, t=20, b=20),
#                             width=FIG_WIDTH+100,
#                             height=FIG_HEIGHT-300,
#                             showlegend=False,
#                             xaxis = dict(showgrid = True,showticklabels = True,titlefont=dict(size=FONT_SIZE),title="Proximity p-value"),
#                             yaxis = dict(showgrid = True,showticklabels = True,titlefont=dict(size=FONT_SIZE),
#                                         title="Embedding score")
#                             )
#     st.plotly_chart(fig_final)
    
                 
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
