import plotly.graph_objects as go
import seaborn as sns
import matplotlib.pyplot as plt


def valid_cells_per_slice(adata, min_genes_list2=None, slice='id'):
    min_genes_list = list(range(0, 800, 100)) if min_genes_list2 is None else min_genes_list2
    library_median_n_genes = adata.obs.groupby(slice)['n_genes'].median().tolist()

    fig = go.Figure()

    s1 = adata.obs.groupby(slice).size()
    s2_list = []
    colors = ['blue', 'red', 'green', 'orange', 'purple', 'pink', 'gray', 'cyan']
    colors_alpha = ['rgba(0, 0, 255, 0.25)', 'rgba(255, 0, 0, 0.25)', 'rgba(0, 255, 0, 0.25)',
                    'rgba(255, 165, 0, 0.25)', 'rgba(128, 0, 128, 0.25)', 'rgba(255, 192, 203, 0.25)',
                    'rgba(128, 128, 128, 0.25)', 'rgba(0, 255, 255, 0.25)']
    for i, min_genes in enumerate(min_genes_list):
        s2 = adata.obs.loc[adata.obs['n_genes'] >= min_genes].groupby(slice).size()
        if s2.empty or not any(s2):
            print(f"Warning: No cells with n_genes >= {min_genes}. Skipping min_genes={min_genes}.")
            continue
        s2_list.append(s2)
        p = s2 / s1
        index_mapping = list(range(len(s1.index)))

        _, ax = plt.subplots()
        sns.regplot(x=library_median_n_genes, y=p, ax=ax, ci=95, scatter=False, color=colors[i % len(colors)])
        x_regplot, y_regplot = ax.lines[0].get_data()
        conf_int = ax.collections[0].get_paths()[0]

        fig.add_trace(go.Scatter(
            x=[library_median_n_genes[index_mapping[j]] for j in range(len(index_mapping))],
            y=p,
            mode='markers',
            marker=dict(size=5, color=colors[i % len(colors)]),
            name=f'min_genes={min_genes}'
        ))

        fig.add_trace(go.Scatter(
            x=x_regplot,
            y=y_regplot,
            mode='lines',
            line=dict(width=1, color=colors[i % len(colors)]),
            name=f'Regression min_genes={min_genes}',
            showlegend=False
        ))

        x_conf_int, y_conf_int = conf_int.vertices.T
        fig.add_trace(go.Scatter(
            x=x_conf_int,
            y=y_conf_int,
            mode='lines',
            line=dict(width=1, color='rgba(68, 122, 219, 0)'),
            name=f'Confidence Interval min_genes={min_genes}',
            fill='tozeroy',
            fillcolor=colors_alpha[i % len(colors_alpha)],
            showlegend=False
        ))

        fig.update_layout(
            title=dict(text='Valid cell ratio post min_genes', x=0.5, y=0.9),
            yaxis_title='Proportion of valid cells',
            xaxis_title='Median ngenes per slice',
            plot_bgcolor='white',
            legend=dict(
                orientation="v",
                x=1.05,
                y=0.5
            )
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


def create_bar_chart(adata, min_genes_list2=None, slice='id'):
    min_genes_list = list(range(0, 800, 100)) if min_genes_list2 is None else min_genes_list2
    fig = go.Figure()

    colors = ['blue', 'red', 'green', 'orange', 'purple', 'pink', 'gray', 'cyan']

    s1 = adata.obs.groupby(slice).size()

    for i, min_genes in enumerate(min_genes_list):
        s2 = adata.obs.loc[adata.obs['n_genes'] >= min_genes].groupby(slice).size()
        if s2.empty or not any(s2):
            print(f"Warning: No cells with n_genes >= {min_genes}. Skipping min_genes={min_genes}.")
            continue
        p = s2 / s1

        fig.add_trace(go.Bar(
            x=list(s1.index),
            y=p,
            name=f'min_genes={min_genes}',
            marker_color=colors[i % len(colors)]
        ))

    fig.update_layout(
        title=dict(text='Valid cell ratio post min_genes', x=0.5, y=0.9),
        barmode='group',
        yaxis_title='Proportion of valid cells',
        xaxis_title='Slice',
        plot_bgcolor='white',
        legend=dict(
            orientation="v",
            x=1.05,
            y=0.5
        )
    )

    return {'data': fig['data'], 'layout': fig['layout']}


def create_plot1(adata, slice='id', min_genes_list2=None):
    data_fig = valid_cells_per_slice(adata, slice=slice, min_genes_list2=min_genes_list2)
    return data_fig


def create_plot2(adata, slice='id', min_genes_list2=None):
    data_fig = create_bar_chart(adata, slice=slice, min_genes_list2=min_genes_list2)
    return data_fig