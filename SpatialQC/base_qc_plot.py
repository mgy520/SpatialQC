import plotly.graph_objects as go
from plotly.subplots import make_subplots
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


def layout(fig):
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
    return fig


def plot_distribution(adata, obs='n_genes'):
    fig = go.Figure()

    data = adata.obs[obs]
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

    fig.add_trace(go.Box(y=data2,
                         q1=[q1],
                         median=[median],
                         q3=[q3],
                         lowerfence=[lowerfence],
                         upperfence=[upperfence],
                         line_color='black',
                         fillcolor='rgba(220, 65, 80, 0.5)',
                         name='',
                         hoverinfo='y',
                         opacity=0.7))

    fig.update_layout(title=dict(text=obs, x=0.5, y=0.9), plot_bgcolor='white', yaxis_title=obs)
    layout(fig)

    return {'data': fig['data'], 'layout': fig['layout']}


def plot_kde(adata, obs='n_genes'):
    tmp = adata.obs[obs]
    kde = sns.kdeplot(tmp, log_scale=True, fill=False)
    x_vals = []
    y_vals = []
    if len(kde.get_lines()) > 0:
        x_vals = kde.get_lines()[0].get_xdata()
        y_vals = kde.get_lines()[0].get_ydata()
    plt.clf()

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(x=x_vals.copy(), y=y_vals.copy(), mode='lines', line=dict(color='blue'), name='Kde')
    )
    max_y = max(y_vals) if len(y_vals) > 0 else 0
    if obs == 'n_counts':
        fig.add_trace(
            go.Scatter(x=[500, 500], y=[0, max_y], mode='lines', line=dict(color='red'),
                       name='n_counts=500')
        )
    elif obs == 'n_genes':
        fig.add_trace(
            go.Scatter(x=[300, 300], y=[0, max_y], mode='lines', line=dict(color='red'),
                       name='n_genes=300')
        )
    else:
        fig.add_trace(
            go.Scatter(x=[0.8, 0.8], y=[0, max_y], mode='lines', line=dict(color='red'),
                       name='log10GenesPerUMI=0.8')
        )
    if obs == 'log10GenesPerUMI':
        fig.update_layout(
            xaxis_title=obs,
            yaxis_title='cell density',
            plot_bgcolor='white',
            legend=dict(
                x=0.55,
                y=1,
                traceorder='normal',
                orientation='v',
                bgcolor='rgba(0,0,0,0)'
            ),
            yaxis=dict(
                range=[0, max_y+3],
            ),
        )
    else:
        fig.update_layout(
            xaxis_title=obs,
            yaxis_title='log10 cell density',
            plot_bgcolor='white',
            legend=dict(
                x=0,
                y=1,
                traceorder='normal',
                orientation='v',
                bgcolor='rgba(0,0,0,0)'
            ),
            xaxis=dict(
                type='log',
                range=[1, 5]
            ),
            yaxis=dict(
                range=[0, max_y + 0.05],
            ),
        )
    layout(fig)

    return {'data': fig['data'], 'layout': fig['layout']}


def mito_ratio(adata, mito_percent=0.1):
    n_counts = adata.obs['n_counts']
    n_genes = adata.obs['n_genes']
    percent_mt = adata.obs['percent_mt']
    fig = go.Figure()

    fig.add_trace(
        go.Scattergl(
            x=n_counts,
            y=n_genes,
            mode='markers',
            marker=dict(
                size=4,
                color=percent_mt,
                colorscale='Greys',
                colorbar=dict(title='Percent mt'),
                cmin=0,
                cmax=mito_percent,
                symbol='circle',
            ),
            text='Ratio:' + percent_mt.astype(str),
            hoverinfo='text',
            name=f'percent_mt<={mito_percent}'
        )
    )

    fig.add_shape(
        type='line',
        x0=0,
        x1=n_counts.max(),
        y0=250,
        y1=250,
    )

    fig.add_shape(
        type='line',
        x0=500,
        x1=500,
        y0=0,
        y1=n_genes.max(),
    )

    fig.add_annotation(
        x=n_counts.max(),
        y=250,
        text='n_genes = 250',
        showarrow=True,
        arrowhead=4,
        ax=-50,
        ay=-30,
    )

    fig.add_annotation(
        x=500,
        y=n_genes.max(),
        text='n_counts = 500',
        showarrow=True,
        arrowhead=4,
        ax=-30,
        ay=-50,
    )
    high_percent_points = adata[adata.obs['percent_mt'] > mito_percent]
    fig.add_trace(
        go.Scattergl(
            x=high_percent_points.obs['n_counts'],
            y=high_percent_points.obs['n_genes'],
            mode='markers',
            marker=dict(
                size=4,
                color='red',
                symbol='circle',
            ),
            text='Percent Mito: ' + high_percent_points.obs['percent_mt'].astype(str),
            hoverinfo='text',
            name=f'percent_mt>{mito_percent}'
        )
    )

    fig.update_layout(
        title=dict(text='n_genes~n_counts colored by percent_mt', x=0.5, y=0.92, font=dict(size=25)),
        xaxis_title=dict(text='n_counts', font=dict(size=20)),
        plot_bgcolor='white',
        yaxis_title=dict(text='n_genes', font=dict(size=20)),
        legend=dict(yanchor="bottom", xanchor="right", x=0.95, y=0.01, font=dict(size=20),
                    bgcolor='rgba(255, 255, 255, 0)'),
    )
    layout(fig)

    return {'data': fig['data'], 'layout': fig['layout']}

def kde_dis(adata2, obs='n_genes'):
    fig2 = plot_distribution(adata=adata2, obs=obs)
    fig1 = plot_kde(adata=adata2, obs=obs)
    fig = make_subplots(rows=1, cols=2)
    for trace in fig1['data']:
        fig.add_trace(trace, row=1, col=1)
    for trace in fig2['data']:
        fig.add_trace(trace, row=1, col=2)

    fig.update_xaxes(fig1['layout']['xaxis'], row=1, col=1)
    fig.update_yaxes(fig1['layout']['yaxis'], row=1, col=1)
    fig.update_xaxes(fig2['layout']['xaxis'], row=1, col=2)
    fig.update_yaxes(fig2['layout']['yaxis'], row=1, col=2)
    fig.update_layout(plot_bgcolor='white')

    return {'data': fig['data'], 'layout': fig['layout']}
