import argparse


def main():
    def str2bool(v):
        if isinstance(v, bool):
            return v
        if v == 'True':
            return True
        elif v == 'False':
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected. Please enter True or False.')


    def mround(match):
        return "{:.4f}".format(float(match.group()))

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    group0 = parser.add_argument_group('Required parameters')
    group4 = parser.add_argument_group('Markers input options (From CellMarker2.0)')
    group5 = parser.add_argument_group('Other parameters related to images in html reports')
    group1 = parser.add_argument_group('Modify the scoring rules for Slice Scores')
    group2 = parser.add_argument_group('Filter related parameters')
    group3 = parser.add_argument_group('Other options')

    group0.add_argument('--adata', metavar='', help='The path to your adata file.')
    group0.add_argument('--markers', metavar='', help='The path to your markers csv file. Prior knowledge,'
                        ' or scRNA-seq of the same tissue to find differential genes, or you may be able to find the '
                                                      'appropriate ones from the Markers input options below')
    group4.add_argument('--species', metavar='', help='The species of your sample.')
    group4.add_argument('--tissue_class', metavar='', help='The tissue class of your sample.')
    group4.add_argument('--tissue_type', metavar='', help='The tissue type of your sample.')
    group4.add_argument('--cancer_type', default='Normal', metavar='', help='The cancer type of your sample.'
                                                                            '\ndefault: Normal')
    group0.add_argument('--slice_number', type=str, default='multiple', metavar='',
                        help='the number of slices. multiple or 1.\ndefault: multiple')
    group0.add_argument('--slice', default='id', metavar='', help="The name that represents the slice "
                        "identifier in adata.obs.If there is only one slice, ignore this parameter.\ndefault: 'id'")
    group0.add_argument('--mito', default='Mt-', metavar='',
                        help="The pattern of mitochondrial genes.\ndefault: 'Mt-'")
    group0.add_argument('--doublet', type=str2bool, default=True, metavar='',
                        help="Whether to identify doublet cells. Recommended for data with single-cell resolution\ndefault: True")
    group2.add_argument('--mito_percent', default='0.1', metavar='', type=float,
                        help="Filter cells with mitochondrial proportion higher than mito_percent.\ndefault: 0.1")
    group0.add_argument('--ribo', default='Rps, Rpl', metavar='',
                        help="The pattern of ribosome genes.\ndefault: 'Rps, Rpl'")
    group0.add_argument('--hemo', default='Hbb, Hba', metavar='',
                        help="The pattern of hemoglobin genes.\ndefault: 'Hbb, Hba'")
    group5.add_argument('--bin_value', default=100, metavar='', type=int,
                        help='Values of n_genes bin intervals applied to the Marker Proportion button.\ndefault: 100')
    group5.add_argument('--min_genes_list', default=list(range(0, 1200, 100)), metavar='', nargs='+', type=int,
                        help='--min_genes_list 100 200 350 ...\n'
                             'used for html buttons:\nCell Number Post Filter\nMarkers '
                             'Proportion Post Filter\nMarkers Detected Post Filter\ndefault: list(range(0, 1200, 100))')
    group5.add_argument('--min_genes_list2', default=list(range(0, 800, 100)), metavar='', nargs='+', type=int,
                        help='--min_genes_list2 100 200 350 ...\nused for button:\nValid Cell Post min_genes.\ndefault: '
                             'list(range(0, 800, 100))')
    group5.add_argument('--min_cells_list', default=[3, 10, 20], metavar='', nargs='+', type=int,
                        help='--min_cells_list 3 10 20 ...\nused for html buttons:\nMarkers Proportion '
                             'Post Filter\nMarkers Detected Post Filter\ndefault: [3, 10, 20]')
    group2.add_argument('--f', type=str2bool, metavar='', default=True, help='Whether to filter adata.\ndefault: True')
    group2.add_argument('--s', default=5, metavar='', type=int,
                        help='Sections with a median score less than s will be removed.\ndefault: 5')
    group2.add_argument('--min_genes', default=None, metavar='',help='Provide your min_genes, otherwise determined by --n.\ndefault: None')
    group2.add_argument('--n', default=0.7, metavar='', type=float,
                        help='Determine the value of min_genes to ensure that the valid cell ratio is greater than n.'
                             '\ndefault: 0.7')
    group2.add_argument('--min_cells', default=None, metavar='',help='Provide your min_cells, otherwise determined by --l.\ndefault: None')
    group2.add_argument('--l', default=0.99, metavar='', type=float,
                        help='After filtering cells, determine the value of min_cells to ensure that the proportion '
                             'of marker genes is greater than l among the remaining detected markers.\ndefault: 0.99')
    group1.add_argument('--s1', type=float, nargs=2, default=[-1, 0], metavar='',
                        help='percent.mt_score (default: -1 0).\npercent_mito of cells > mito_percent: -1\n'
                             'percent_mito of cells <= mito_percent: 0')
    group1.add_argument('--s2', type=float, nargs=3, default=[0.8, 0, 1], metavar='',
                        help='log10GenesPerUMI_score (default: 0.8, 0, 1).\nlog10GenesPerUMI of cells < 0.8: 0\n'
                             'log10GenesPerUMI of cells >= 0.8: 1')
    group1.add_argument('--s3', type=float, nargs=5, default=[0.2, 0.5, 0, 1, 2], metavar='',
                        help='n_genes_score (default: 0.2 0.5 0 1 2).\nn_genes ranking below 20th percentile: 0\n'
                             'between 20th and 50th percentile: 1\nabove 50th percentile: 2')
    group1.add_argument('--s4', type=float, nargs=7, default=[0.2, 0.5, 0.8, 0, 1, 2, 3], metavar='',
                        help='markerDetectionRatio_score (default: 0.2, 0.5, 0.8, 0, 1, 2, 3).\n'
                             'markerDetectionRatio ranking below 20th percentile: 0\n'
                             'between 20th and 50th percentile: 1\nbetween 50th and 80th percentile: 2\nabove 80th percentile: 3')
    group1.add_argument('--s5', type=float, nargs=7, default=[0.2, 0.5, 0.8, 0, 1, 2, 3], metavar='',
                        help='markerProportion_score (default: 0.2, 0.5, 0.8, 0, 1, 2, 3).\n'
                             'markerProportion_score ranking below 20th percentile: 0\n'
                             'between 20th and 50th percentile: 1\nbetween 50th and 80th percentile: 2\nabove 80th percentile: 3')
    group1.add_argument('--s6', type=float, nargs=7, default=[0.2, 0.5, 0.8, 0, 1, 2, 3], metavar='',
                        help='markerCountsRatio_score (default: 0.2, 0.5, 0.8, 0, 1, 2, 3).\n'
                             'markerCountsRatio_score ranking below 20th percentile: 0\n'
                             'between 20th and 50th percentile: 1\nbetween 50th and 80th percentile: 2\nabove 80th percentile: 3')
    group1.add_argument('--s7', type=float, nargs=2, default=[-4, 0], metavar='',
                        help='doublet_score (default: -4 0).\ndoublet cells: -4\nnot doublet cells: 0')
    group1.add_argument('--s8', type=float, nargs=5, default=[0.2, 0.5, 0, 1, 2], metavar='',
                        help='n_counts_score (default: 0.2 0.5 0 1 2).\nn_counts ranking below 20th percentile: 0\n'
                             'between 20th and 50th percentile: 1\nabove 50th percentile: 2')
    group3.add_argument('--output', default='./', metavar='', help="output directory.\ndefault: './'")
    group3.add_argument('--o1', type=str, metavar='', default='report.html',
                        help='the filename for the output of html report.\ndefault: report.html')
    group3.add_argument('--o2', type=str, metavar='', default='filtered.h5ad',
                        help='the filename for the output of filtered adata.\ndefault: filtered.h5ad)')
    group3.add_argument('--j', type=int, metavar='', default=8,
                        help='The maximum number of concurrently running jobs. If set to 1, parallelism is not used. If set to -1,'
                             ' all CPUs are used. For n_jobs less than -1, (n_cpus + 1 + n_jobs) CPUs are used.\ndefault: 8')
    args = parser.parse_args()

    import scanpy as sc
    from . import basic_statistics
    from . import slice_score
    from . import markers_detected_by_ngenes_cellbin
    from . import markers_detected_by_filter_cells_genes
    from . import markers_detected_by_filter_cells_genes2
    from . import valid_cells
    from . import base_qc_plot


    adata = sc.read_h5ad(args.adata)
    sc.pp.filter_cells(adata, min_genes=1)
    sc.pp.filter_cells(adata, min_counts=1)
    sc.pp.filter_genes(adata, min_cells=1)
    ribo_genes = adata.var_names.str.startswith(tuple(args.ribo.split(',')))
    hb_genes = adata.var_names.str.startswith(tuple(args.hemo.split(',')))
    import numpy as np
    adata.obs['percent_hb'] = np.sum(
        adata[:, hb_genes].X, axis=1) / np.sum(adata.X, axis=1)
    adata.obs['percent_ribo'] = np.sum(
        adata[:, ribo_genes].X, axis=1) / np.sum(adata.X, axis=1)
    markers = None
    import pandas as pd
    if args.markers is not None:
        markers = pd.read_csv(args.markers).squeeze("columns").tolist()

    button_plots = []
    slice_score_fig_all, modified_adata = slice_score.create_plot(
        adata, markers, args.species, args.tissue_class, args.tissue_type, args.cancer_type, slice=args.slice,
        mito=args.mito, mito_p=args.mito_percent, doublet=args.doublet, s1=args.s1, s2=args.s2, s3=args.s3, s4=args.s4, s5=args.s5,
        s6=args.s6, s7=args.s7, s8=args.s8)
    del adata
    if modified_adata.obs[args.slice].nunique() > 6:
        valid_cell_fig_all = valid_cells.create_plot1(modified_adata, slice=args.slice, min_genes_list2=args.min_genes_list2)
    else:
        valid_cell_fig_all = valid_cells.create_plot2(modified_adata, slice=args.slice, min_genes_list2=args.min_genes_list2)
    from joblib import Parallel, delayed
    functions_and_params = [
        (basic_statistics.create_plot,
         (modified_adata, markers, args.species, args.tissue_class, args.tissue_type, args.cancer_type, args.slice)),
        (base_qc_plot.kde_dis, (modified_adata, 'n_counts')),
        (base_qc_plot.kde_dis, (modified_adata, 'n_genes')),
        (base_qc_plot.kde_dis, (modified_adata, 'log10GenesPerUMI')),
        (base_qc_plot.plot_distribution, (modified_adata, 'percent_mt')),
        (base_qc_plot.mito_ratio, (modified_adata, args.mito_percent)),
        (base_qc_plot.plot_distribution, (modified_adata, 'percent_ribo')),
        (base_qc_plot.plot_distribution, (modified_adata, 'percent_hb')),
        (markers_detected_by_ngenes_cellbin.create_plot, (
        modified_adata, markers, args.species, args.tissue_class, args.tissue_type, args.cancer_type, args.bin_value)),
        (basic_statistics.cell_numbers_by_min_genes, (modified_adata, args.min_genes_list)),
        (markers_detected_by_filter_cells_genes2.create_plot, (
        modified_adata, markers, args.species, args.tissue_class, args.tissue_type, args.cancer_type,
        args.min_genes_list, args.min_cells_list)),
        (markers_detected_by_filter_cells_genes.create_plot, (
        modified_adata, markers, args.species, args.tissue_class, args.tissue_type, args.cancer_type,
        args.min_genes_list, args.min_cells_list))
    ]
    results = Parallel(n_jobs=args.j)(delayed(func)(*params) for func, params in functions_and_params)
    (basic_statistics_fig_all, n_counts_per_cell_fig_all, n_genes_per_cell_fig_all, log10_genes_per_umi_fig_all, mito_box_fig_all,
     mito_ratio_fig_all, ribo_fig_all, hb_fig_all, marker_proportion_fig_all, cell_number_fig_all, marker_post_min_genes_fig_all,
     marker_detected_post_filter_fig_all) = results

    button_plots.append({
        'basic_statistics_all': {
            'data': basic_statistics_fig_all['data'],
            'layout': basic_statistics_fig_all['layout']
        },
        'slice_score_all': {
            'data': slice_score_fig_all['data'],
            'layout': slice_score_fig_all['layout']
        },
        'n_counts_per_cell_all': {
            'data': n_counts_per_cell_fig_all['data'],
            'layout': n_counts_per_cell_fig_all['layout']
        },
        'n_genes_per_cell_all': {
            'data': n_genes_per_cell_fig_all['data'],
            'layout': n_genes_per_cell_fig_all['layout']
        },
        'log10GenesPerUMI_all': {
            'data': log10_genes_per_umi_fig_all['data'],
            'layout': log10_genes_per_umi_fig_all['layout']
        },
        'mito_box_all': {
            'data': mito_box_fig_all['data'],
            'layout': mito_box_fig_all['layout']
        },
        'mito_ratio_all': {
            'data': mito_ratio_fig_all['data'],
            'layout': mito_ratio_fig_all['layout']
        },
        'ribo_all': {
            'data': ribo_fig_all['data'],
            'layout': ribo_fig_all['layout']
        },
        'hb_all': {
            'data': hb_fig_all['data'],
            'layout': hb_fig_all['layout']
        },
        'marker_proportion_all': {
            'data': marker_proportion_fig_all['data'],
            'layout': marker_proportion_fig_all['layout']
        },
        'valid_cell_all': {
            'data': valid_cell_fig_all['data'],
            'layout': valid_cell_fig_all['layout']
        },
        'cell_number_all': {
            'data': cell_number_fig_all['data'],
            'layout': cell_number_fig_all['layout']
        },
        'marker_post_min_cells_all': {
            'data': marker_post_min_genes_fig_all['data'],
            'layout': marker_post_min_genes_fig_all['layout']
        },
        'marker_detected_post_filter_all': {
            'data': marker_detected_post_filter_fig_all['data'],
            'layout': marker_detected_post_filter_fig_all['layout']
        }
    })

    # for slice
    def process_slice(slice_id):
        bdata = modified_adata[modified_adata.obs[args.slice] == slice_id].copy()
        sc.pp.filter_genes(bdata, min_cells=1)
        basic_statistics_fig_slice = basic_statistics.create_plot(
            bdata, markers, args.species, args.tissue_class, args.tissue_type, args.cancer_type, slice=args.slice)
        s1_fig_slice = slice_score.specific_score(bdata, 'doublet_score')
        s2_fig_slice = slice_score.specific_score(bdata, 'percent.mt_score')
        s3_fig_slice = slice_score.specific_score(bdata, 'log10GenesPerUMI_score')
        s4_fig_slice = slice_score.specific_score(bdata, 'n_genes_score')
        s5_fig_slice = slice_score.specific_score(bdata, 'markerDetectionRatio_score')
        s6_fig_slice = slice_score.specific_score(bdata, 'markerProportion_score')
        s7_fig_slice = slice_score.specific_score(bdata, 'markerCountsRatio_score')
        s8_fig_slice = slice_score.specific_score(bdata, 'n_counts_score')
        n_counts_per_cell_fig_slice = base_qc_plot.kde_dis(bdata, obs='n_counts')
        n_genes_per_cell_fig_slice = base_qc_plot.kde_dis(bdata, obs='n_genes')
        log10_genes_per_umi_fig_slice = base_qc_plot.kde_dis(bdata, obs='log10GenesPerUMI')
        mito_box_fig_slice = base_qc_plot.plot_distribution(bdata, obs='percent_mt')
        mito_ratio_fig_slice = base_qc_plot.mito_ratio(bdata, mito_percent=args.mito_percent)
        ribo_fig_slice = base_qc_plot.plot_distribution(bdata, obs='percent_ribo')
        hb_fig_slice = base_qc_plot.plot_distribution(bdata, obs='percent_hb')
        marker_proportion_fig_slice = markers_detected_by_ngenes_cellbin.create_plot(
            bdata, markers, args.species, args.tissue_class, args.tissue_type,
            args.cancer_type, bin_value=args.bin_value)
        cell_number_fig_slice = basic_statistics.cell_numbers_by_min_genes(bdata, min_genes_list=args.min_genes_list)
        valid_cell_fig_slice = valid_cells.create_plot2(bdata, slice=args.slice)
        marker_post_min_genes_fig_slice = markers_detected_by_filter_cells_genes2.create_plot(
            bdata, markers, args.species, args.tissue_class, args.tissue_type, args.cancer_type,
            min_genes_list=args.min_genes_list, min_cells_list=args.min_cells_list)
        marker_detected_post_filter_fig_slice = markers_detected_by_filter_cells_genes.create_plot(
            bdata, markers, args.species, args.tissue_class, args.tissue_type, args.cancer_type,
            min_genes_list=args.min_genes_list, min_cells_list=args.min_cells_list)

        return {
            'basic_statistics_' + slice_id: {
                'data': basic_statistics_fig_slice['data'],
                'layout': basic_statistics_fig_slice['layout']
            },
            's1_score_' + slice_id: {
                'data': s1_fig_slice['data'],
                'layout': s1_fig_slice['layout']
            },
            's2_score_' + slice_id: {
                'data': s2_fig_slice['data'],
                'layout': s2_fig_slice['layout']
            },
            's3_score_' + slice_id: {
                'data': s3_fig_slice['data'],
                'layout': s3_fig_slice['layout']
            },
            's4_score_' + slice_id: {
                'data': s4_fig_slice['data'],
                'layout': s4_fig_slice['layout']
            },
            's5_score_' + slice_id: {
                'data': s5_fig_slice['data'],
                'layout': s5_fig_slice['layout']
            },
            's6_score_' + slice_id: {
                'data': s6_fig_slice['data'],
                'layout': s6_fig_slice['layout']
            },
            's7_score_' + slice_id: {
                'data': s7_fig_slice['data'],
                'layout': s7_fig_slice['layout']
            },
            's8_score_' + slice_id: {
                'data': s8_fig_slice['data'],
                'layout': s8_fig_slice['layout']
            },
            'n_counts_per_cell_' + slice_id: {
                'data': n_counts_per_cell_fig_slice['data'],
                'layout': n_counts_per_cell_fig_slice['layout']
            },
            'n_genes_per_cell_' + slice_id: {
                'data': n_genes_per_cell_fig_slice['data'],
                'layout': n_genes_per_cell_fig_slice['layout']
            },
            'log10GenesPerUMI_' + slice_id: {
                'data': log10_genes_per_umi_fig_slice['data'],
                'layout': log10_genes_per_umi_fig_slice['layout']
            },
            'mito_box_' + slice_id: {
                'data': mito_box_fig_slice['data'],
                'layout': mito_box_fig_slice['layout']
            },
            'mito_ratio_' + slice_id: {
                'data': mito_ratio_fig_slice['data'],
                'layout': mito_ratio_fig_slice['layout']
            },
            'ribo_' + slice_id: {
                'data': ribo_fig_slice['data'],
                'layout': ribo_fig_slice['layout']
            },
            'hb_' + slice_id: {
                'data': hb_fig_slice['data'],
                'layout': hb_fig_slice['layout']
            },
            'marker_proportion_' + slice_id: {
                'data': marker_proportion_fig_slice['data'],
                'layout': marker_proportion_fig_slice['layout']
            },
            'valid_cell_' + slice_id: {
                'data': valid_cell_fig_slice['data'],
                'layout': valid_cell_fig_slice['layout']
            },
            'cell_number_' + slice_id: {
                'data': cell_number_fig_slice['data'],
                'layout': cell_number_fig_slice['layout']
            },
            'marker_post_min_cells_' + slice_id: {
                'data': marker_post_min_genes_fig_slice['data'],
                'layout': marker_post_min_genes_fig_slice['layout']
            },
            'marker_detected_post_filter_' + slice_id: {
                'data': marker_detected_post_filter_fig_slice['data'],
                'layout': marker_detected_post_filter_fig_slice['layout']
            }
        }
    if args.slice_number == 'multiple':
        slice_ids = modified_adata.obs[args.slice].unique()
        results = Parallel(n_jobs=args.j)(
            delayed(process_slice)(slice_id) for slice_id in modified_adata.obs[args.slice].unique())
        button_plots.extend(results)
        slice_names = list(slice_ids)

    # report
    import json
    import os
    import plotly
    import re
    import gzip
    import base64
    js_data = json.dumps(button_plots, cls=plotly.utils.PlotlyJSONEncoder, separators=(',', ':'))
    js_data = re.sub(re.compile(r"\d+\.\d{8,}"), mround, js_data)
    js_data2 = gzip.compress(js_data.encode('utf-8'))
    js_data3 = base64.b64encode(js_data2).decode('utf-8')

    slice_names_json = json.dumps(['all'] + slice_names)
    dir1 = os.path.dirname(os.path.realpath(__file__))
    html_path = os.path.join(dir1, 'templates', 'report_template.html')

    with open(html_path, 'r') as f:
        template = f.read()
    html = template.replace('{{data}}', js_data3)
    html = html.replace('{{slice_names}}', f'{slice_names_json}')

    os.makedirs(args.output, exist_ok=True)
    output_directory = args.output
    with open(os.path.join(output_directory, args.o1), 'w') as f:
        f.write(html)

    # filter
    if args.f:
        selected_slices = modified_adata.obs.groupby(args.slice)['cell_score'].median() > args.s
        cdata = modified_adata[modified_adata.obs[args.slice].isin(selected_slices[selected_slices].index)].copy()
        if cdata.obs['predicted_doublets'].isnull().any():
            print("Warning: 'predicted_doublets' contains null values. Maybe you need to double-check doublet")
        cdata = cdata[cdata.obs['predicted_doublets'].isin([False, None])].copy()
        if args.min_genes is None or args.min_genes == "None":
            thresholds = cdata.obs.groupby(args.slice).apply(lambda group: group['n_genes'].quantile(1 - args.n))
            final_threshold = int(thresholds.min() / 10) * 10
            print(f"Automatic threshold for n_genes is: {final_threshold}")
        else:
            final_threshold = int(args.min_genes)
            print(f"Using user-provided min_genes: {final_threshold}")
        sc.pp.filter_cells(cdata, min_genes=final_threshold)
        cdata = cdata[cdata.obs['percent_mt'] <= args.mito_percent].copy()
        sc.pp.filter_genes(cdata, min_cells=1)
        if args.min_cells is None or args.min_cells == "None":
            from .get_markers import get_markers_from_db
            markers = get_markers_from_db(args.species, args.tissue_class, args.tissue_type,
                                          args.cancer_type) if args.markers is None else pd.read_csv(
                args.markers).squeeze(
                "columns").tolist()
            i = list(set(markers).intersection(cdata.var_names))
            selected = cdata.var.loc[i]
            threshold = int(selected['n_cells'].quantile(1 - args.l))
            print(f"The suitable threshold for n_cells is: {threshold}")
        else:
            threshold = int(args.min_cells)
            print(f"Using user-provided min_cells: {threshold}")
        sc.pp.filter_genes(cdata, min_cells=threshold)
        cdata.obs = cdata.obs.drop(['percent.mt_score', 'log10GenesPerUMI_score', 'n_genes_score', 'markerDetectionRatio_score',
         'markerProportion_score', 'markerCountsRatio_score', 'doublet_score', 'n_counts_score', 'cell_score', 'remain_ratio',
                                    'marker_ratio', 'marker_exp', 'predicted_doublets'], axis=1)
        cdata.write(filename=os.path.join(output_directory, args.o2))


if __name__ == '__main__':
    main()
