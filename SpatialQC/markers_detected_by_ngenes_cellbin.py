import plotly.graph_objects as go
import numpy as np
import scipy
from .get_markers import marker_decorator


@marker_decorator
def markers_detected_by_diff_ngenes_bar(adata, markers, bin_value=100):
    bins = np.arange(0, np.ceil(adata.obs['n_genes'].max() / bin_value) * bin_value, bin_value).astype(int)
    adata.obs['n_genes_subset'] = np.digitize(adata.obs['n_genes'], bins)
    markers_set = set(markers)

    counts = []
    for i in range(len(adata)):
        nonzero_indices = adata.X[i].nonzero()
        if isinstance(adata.X, scipy.sparse.csr_matrix):
            genes = adata.var_names[nonzero_indices[1]]
        else:
            genes = adata.var_names[nonzero_indices]
        intersect = set(genes).intersection(markers_set)
        counts.append(len(intersect))

    adata.obs["Markers intersect counts"] = counts
    subset_counts = (adata.obs.groupby('n_genes_subset')['Markers intersect counts'].mean() / len(markers_set)).tolist()

    fig = go.Figure(data=go.Bar(x=[f"{bins[i]}-{bins[i+1]}" for i in range(len(bins) - 1)], y=subset_counts))
    fig.update_layout(title=dict(text='Marker proportion ~ n_genes',x=0.5, y=0.9),
                      xaxis_title='n_genes',
                      yaxis_title='Markers proportion',
                      plot_bgcolor='white')
    fig.update_xaxes(
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
        gridcolor='lightgrey'
    )
    fig.update_yaxes(
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
        gridcolor='lightgrey'
    )

    return {'data': fig['data'], 'layout': fig['layout']}


def create_plot(adata, markers, species=None, tissue_class=None, tissue_type=None, cancer_type=None, bin_value=100):
    data_fig = markers_detected_by_diff_ngenes_bar(adata, markers, species, tissue_class, tissue_type, cancer_type, bin_value=bin_value)
    return data_fig
