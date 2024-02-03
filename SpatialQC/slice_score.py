import pandas as pd
import numpy as np
import scrublet as scr
import plotly.graph_objects as go
from .get_markers import marker_decorator
from scipy.sparse import issparse


@marker_decorator
def slice_score(adata, markers, mito='Mt-', slice='id', mito_p=0.1, doublet=True, s1=None, s2=None,
                s3=None, s4=None, s5=None, s6=None, s7=None, s8=None):
    default_values = {
        's1': [-1, 0],
        's2': [0.8, 0, 1],
        's3': [0.2, 0.5, 0, 1, 2],
        's4': [0.2, 0.5, 0.8, 0, 1, 2, 3],
        's5': [0.2, 0.5, 0.8, 0, 1, 2, 3],
        's6': [0.2, 0.5, 0.8, 0, 1, 2, 3],
        's7': [-4, 0],
        's8': [0.2, 0.5, 0, 1, 2]
    }
    for var_name, default_value in default_values.items():
        if locals()[var_name] is None:
            locals()[var_name] = default_value

    mito_genes = adata.var_names.str.startswith(mito)

    adata.obs['percent_mt'] = np.sum(adata[:, mito_genes].X, axis=1) / np.sum(adata.X, axis=1)
    adata.obs['log10GenesPerUMI'] = np.log10(adata.obs['n_genes']) / np.log10(adata.obs['n_counts'])

    genes_b = set(markers)
    adata_genes = adata.var_names
    common_genes = adata_genes.isin(genes_b)
    adata.obs['remain_ratio'] = np.sum(adata[:, common_genes].X > 0, axis=1) / len(genes_b)
    adata.obs['marker_ratio'] = np.sum(adata[:, common_genes].X > 0, axis=1) / np.sum(adata.X > 0, axis=1)

    if issparse(adata.X):
        adata.obs['marker_exp'] = np.sum(adata[:, common_genes].X, axis=1).A1 / adata.obs['n_counts']
    else:
        adata.obs['marker_exp'] = np.sum(adata[:, common_genes].X, axis=1) / adata.obs['n_counts']

    if doublet:
        slice_list = adata.obs[slice].cat.categories.tolist()
        alldata = {}
        for i in slice_list:
            adata_filtered = adata[adata.obs[slice] == i,]
            scrub = scr.Scrublet(adata_filtered.X)
            out = scrub.scrub_doublets(verbose=False)
            alldata[i] = pd.DataFrame({'doublet_score': out[0], 'predicted_doublets': out[1]},
                                      index=adata_filtered.obs.index)

        scr_result = pd.concat(alldata.values())
        adata.obs['scrublet_doublet_scores'] = scr_result['doublet_score']
        adata.obs['predicted_doublets'] = scr_result['predicted_doublets']
    else:
        adata.obs['predicted_doublets'] = False

    adata.obs['percent.mt_score'] = np.where(adata.obs['percent_mt'] > mito_p, s1[0], s1[1])
    adata.obs['log10GenesPerUMI_score'] = np.where(adata.obs['log10GenesPerUMI'] < s2[0], s2[1], s2[2])
    percentiles = adata.obs['n_genes'].quantile([s3[0], s3[1]])
    adata.obs['n_genes_score'] = np.where(adata.obs['n_genes'] <= percentiles[s3[0]], s3[2],
                                          np.where(adata.obs['n_genes'] <= percentiles[s3[1]], s3[3], s3[4]))
    percentiles = adata.obs['remain_ratio'].quantile([s4[0], s4[1], s4[2]])
    adata.obs['markerDetectionRatio_score'] = np.where(adata.obs['remain_ratio'] <= percentiles[s4[0]], s4[3],
                                                       np.where(adata.obs['remain_ratio'] <= percentiles[s4[1]], s4[4],
                                                                np.where(adata.obs['remain_ratio'] <= percentiles[s4[2]],
                                                                         s4[5], s4[6])))
    percentiles = adata.obs['marker_ratio'].quantile([s5[0], s5[1], s5[2]])
    adata.obs['markerProportion_score'] = np.where(adata.obs['marker_ratio'] <= percentiles[s5[0]], s5[3],
                                                   np.where(adata.obs['marker_ratio'] <= percentiles[s5[1]], s5[4],
                                                            np.where(adata.obs['marker_ratio'] <= percentiles[s5[2]], s5[5],
                                                                     s5[6])))
    percentiles = adata.obs['marker_exp'].quantile([s6[0], s6[1], s6[2]])
    adata.obs['markerCountsRatio_score'] = np.where(adata.obs['marker_exp'] <= percentiles[s6[0]], s6[3],
                                                    np.where(adata.obs['marker_exp'] <= percentiles[s6[1]], s6[4],
                                                             np.where(adata.obs['marker_exp'] <= percentiles[s6[2]], s6[5],
                                                                      s6[6])))
    adata.obs['doublet_score'] = np.where(adata.obs['predicted_doublets'], s7[0], s7[1])
    percentiles = adata.obs['n_counts'].quantile([s8[0], s8[1]])
    adata.obs['n_counts_score'] = np.where(adata.obs['n_counts'] <= percentiles[s8[0]], s8[2],
                                           np.where(adata.obs['n_counts'] <= percentiles[s8[1]], s8[3], s8[4]))

    adata.obs['cell_score'] = adata.obs[
        ['percent.mt_score', 'log10GenesPerUMI_score', 'n_genes_score', 'markerDetectionRatio_score',
         'markerProportion_score', 'markerCountsRatio_score', 'doublet_score', 'n_counts_score']].sum(axis=1)

    df = adata.obs[[slice, 'cell_score']]
    fig = go.Figure()
    for slice_id in df[slice].unique():
        data = df['cell_score'][df[slice] == slice_id]
        q1 = np.percentile(data, 25)
        median = np.median(data)
        q3 = np.percentile(data, 75)
        iqr = q3 - q1
        lower_bound = q1 - 1.5 * iqr
        upper_bound = q3 + 1.5 * iqr
        lowerfence = max(min(data), lower_bound)
        upperfence = min(max(data), upper_bound)
        outliers = data[(data < lower_bound) | (data > upper_bound)]

        data2 = [outliers] if outliers.any() else None

        fig.add_trace(go.Box(
            x=[slice_id],
            y=data2,
            q1=[q1],
            median=[median],
            q3=[q3],
            lowerfence=[lowerfence],
            upperfence=[upperfence],
            boxpoints='outliers',
            line_color='black',
            fillcolor='rgba(220, 65, 80, 0.5)',
            hoverinfo='y',
            opacity=0.7,
            name=slice_id
        ))

    fig.update_layout(
        title=dict(text='Slice score distribution', x=0.5, y=0.9),
        yaxis_title='Cell_score',
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

    data_fig = {'data': fig['data'], 'layout': fig['layout']}
    return data_fig, adata


def create_plot(adata, markers, species=None, tissue_class=None, tissue_type=None, cancer_type=None,
                slice='id', mito='Mt-', mito_p=0.1, doublet=True, s1=None, s2=None, s3=None,
                s4=None, s5=None, s6=None, s7=None, s8=None):
    data_fig, modified_adata = slice_score(adata, markers, species, tissue_class,
                                           tissue_type, cancer_type,
                                           slice=slice, mito=mito,
                                           mito_p=mito_p, doublet=doublet, s1=s1, s2=s2, s3=s3, s4=s4,
                                           s5=s5, s6=s6, s7=s7, s8=s8)
    return data_fig, modified_adata


def specific_score(bdata, s_column):
    coords = bdata.obsm['spatial']
    colors = bdata.obs[s_column]

    fig = go.Figure(data=go.Scatter(
        x=coords[:, 0],
        y=coords[:, 1],
        mode='markers',
        marker=dict(
            size=4,
            color=colors,
            colorscale='Viridis',
            colorbar=dict(title=s_column, xanchor="left", x=0.8, thickness=10)
        )
    ))

    fig.update_layout(
        title=dict(text=s_column, x=0.5, y=0.99),
        autosize=False,
        plot_bgcolor='white',
        xaxis=dict(
            scaleanchor="y",
            scaleratio=1,
            constrain="domain"
        ),
        yaxis=dict(
            constrain="domain"
        ),
        margin=dict(
            l=20,
            r=20,
            b=20,
            t=20,
        )
    )
    fig.update_xaxes(
        mirror=True,
        showticklabels=False,
        showline=True,
        linecolor='black',
        gridcolor='lightgrey'
    )
    fig.update_yaxes(
        mirror=True,
        showticklabels=False,
        showline=True,
        linecolor='black',
        gridcolor='lightgrey'
    )
    return {'data': fig['data'], 'layout': fig['layout']}
