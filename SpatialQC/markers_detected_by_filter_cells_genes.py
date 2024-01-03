import plotly.graph_objects as go
import numpy as np
from scipy.sparse import issparse
from .get_markers import marker_decorator


@marker_decorator
def detect_markers_by_filter_cells(adata, markers, min_genes_list=None, min_cells_list=None):
    if min_cells_list is None:
        min_cells_list = [1, 5, 10]
    if min_genes_list is None:
        min_genes_list = list(range(0, 1200, 100))
    markers_set = set(markers)
    marker_overlap = {frac: [] for frac in min_cells_list}
    genes = adata.var_names.tolist()
    counts = len(set(genes).intersection(set(markers)))
    for frac in min_cells_list:
        for min_genes in min_genes_list:
            valid_cells = adata.obs['n_genes'] >= min_genes
            remaining_cells = adata.obs_names[valid_cells].tolist()
            if issparse(adata.X):
                new_n_cells = adata.X[adata.obs_names.isin(remaining_cells),].getnnz(axis=0)
            else:
                new_n_cells = np.count_nonzero(adata.X[adata.obs_names.isin(remaining_cells),], axis=0)
            valid_genes = new_n_cells >= frac
            remaining_genes = adata.var_names[valid_genes].tolist()
            overlap_genes = set(remaining_genes).intersection(markers_set)
            marker_overlap_percent = len(overlap_genes) / counts
            marker_overlap[frac].append(marker_overlap_percent)

    fig = go.Figure()
    
    for frac in min_cells_list:
        fig.add_trace(go.Scatter(x=min_genes_list, y=marker_overlap[frac], mode='lines', name=f'min_cells={frac}'))

    fig.update_layout(
        title=dict(text='Remaining markers post filter', x=0.5, y=0.9),
        legend=dict(yanchor="bottom", xanchor="left", x=0.01, y=0.01, bgcolor='rgba(0,0,0,0)'),
        xaxis_title='min_genes',
        yaxis_title='Percentage',
        plot_bgcolor='white'
    )
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


def create_plot(adata, markers, species=None, tissue_class=None, tissue_type=None,
                cancer_type=None, min_genes_list=None, min_cells_list=None):
    data_fig = detect_markers_by_filter_cells(adata, markers, species, tissue_class, tissue_type, cancer_type,
                                              min_genes_list=min_genes_list, min_cells_list=min_cells_list)
    return data_fig
