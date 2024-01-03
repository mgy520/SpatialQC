from .get_markers import marker_decorator
import plotly.figure_factory as ff
import plotly.graph_objects as go


@marker_decorator
def basic_statistics(adata, markers, slice='id'):
    cell_count = adata.shape[0]
    gene_count = adata.shape[1]
    slice_count = adata.obs[slice].nunique()
    median_n_genes = adata.obs['n_genes'].median()
    mean_n_genes = adata.obs['n_genes'].mean()
    median_n_UMIs = adata.obs['n_counts'].median()
    mean_n_UMIs = adata.obs['n_counts'].mean()
    provided_marker_counts = len(set(markers))
    detected_marker_counts = len(set(markers).intersection(adata.var_names))
    marker_detected_ratio = detected_marker_counts / len(set(markers))
    marker_proportion = detected_marker_counts / gene_count
    doublet_count = adata.obs['predicted_doublets'].sum()

    data = [['Measure', 'Value'],
            ['Slice counts', slice_count],
            ['Cell counts', cell_count],
            ['Gene counts', gene_count],
            ['Median n_genes', median_n_genes],
            ['Mean n_genes', mean_n_genes],
            ['Median UMIs', median_n_UMIs],
            ['Mean UMIs', mean_n_UMIs],
            ['Doublet counts', doublet_count],
            ['Provided marker counts', provided_marker_counts],
            ['Detected marker counts', detected_marker_counts],
            ['Marker detected ratio', marker_detected_ratio],
            ['Marker proportion', marker_proportion]]

    fig = ff.create_table(data, index=True)

    for i in range(len(fig.layout.annotations)):
        fig.layout.annotations[i].font.size = 18

    return {'data': fig['data'], 'layout': fig['layout']}


def create_plot(adata, markers, species=None, tissue_class=None, tissue_type=None, cancer_type=None, slice='id'):
    data_fig = basic_statistics(adata, markers, species, tissue_class, tissue_type, cancer_type, slice=slice)
    return data_fig


def cell_numbers_by_min_genes(adata, min_genes_list=None):
    min_genes_list = min_genes_list if min_genes_list is not None else list(range(0, 1200, 100))
    num_cells = []
    for min_gene in min_genes_list:
        num_cells.append(sum(adata.obs['n_genes'] >= min_gene))

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=min_genes_list, y=num_cells, mode='lines'))

    fig.update_layout(
        title=dict(text='Remaining cell number post min_genes', x=0.5, y=0.9),
        xaxis_title='min_genes',
        yaxis_title='Cell Number',
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
