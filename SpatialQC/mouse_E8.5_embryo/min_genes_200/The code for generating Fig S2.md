### The code for generating Fig S2

```python
#import the required packages 
import scanpy as sc
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

#load data
bdata=sc.read('E8.5-Atlas_run_min200.h5ad')
sc.pp.filter_genes(bdata,min_cells=1)
sc.pp.filter_cells(bdata,min_genes=0)

#Fig S2.A
plt.figure()
tmp=bdata.obs['n_genes']
sns.kdeplot(tmp, fill=False,log_scale=True,color='blue')
plt.axvline(x=300, color='red', linestyle='-')
plt.xlabel('n_genes', fontfamily='Arial', fontsize=20)
plt.ylabel('log10 cell density', fontfamily='Arial', fontsize=20)
plt.xticks(fontsize=12, fontfamily='Arial')
plt.yticks(fontsize=12, fontfamily='Arial')
plt.text(350, 0.8, 'n_genes=300', color='red', fontsize=10, fontfamily='Arial')
plt.savefig('n_genes_density_plot_min200.svg', format='svg',bbox_inches='tight')

#Fig S2.B
sc.pp.filter_cells(bdata,min_counts=0)
plt.figure()
tmp=bdata.obs['n_counts']
sns.kdeplot(tmp, fill=False,log_scale=True,color='blue')
plt.axvline(x=500, color='red', linestyle='-')
plt.xlabel('n_counts', fontfamily='Arial', fontsize=20)
plt.ylabel('log10 cell density', fontfamily='Arial', fontsize=20)
plt.xticks(fontsize=12, fontfamily='Arial')
plt.yticks(fontsize=12, fontfamily='Arial')
plt.text(550, 0.55, 'n_counts=500', color='red', fontsize=10, fontfamily='Arial')
plt.savefig('n_counts_density_plot_min200.svg', format='svg',bbox_inches='tight')

#Fig S2.E
mito_genes = bdata.var_names.str.startswith('mt')
bdata.obs['percent_mt'] = np.sum(bdata[:, mito_genes].X, axis=1) / np.sum(bdata.X, axis=1)
plt.figure()
tmp = bdata.obs['percent_mt']
sns.violinplot(y=tmp,facecolor=(220/255, 65/255, 80/255, 0.5))
plt.title('Violin Plot of Mitochondrial ratio',fontfamily='Arial', fontsize=20)
plt.ylabel('Mitochondrial ratio', fontfamily='Arial', fontsize=20)
plt.xticks(fontsize=12, fontfamily='Arial')
plt.yticks(fontsize=12, fontfamily='Arial')
plt.savefig('mt_ratio_box_plot.svg', format='svg',bbox_inches='tight')

#load result .h5ad (The cell score was preserved in the filtered.h5ad)
adata=sc.read('/rad/mgy/results/path/SpatialQC/tmp.h5ad')

#Fig S2.C
plt.figure(figsize=(12, 3))
sns.boxplot(data=adata.obs, x='slice', y='cell_score', color=(220/255, 65/255, 80/255, 0.5))
plt.title('Slices total score distribution',fontfamily='Arial', fontsize=20)
plt.xlabel('')
plt.ylabel('Cell Score',fontfamily='Arial', fontsize=20)
plt.xticks(rotation=45,fontsize=12, fontfamily='Arial')
plt.yticks(fontsize=12, fontfamily='Arial')
plt.savefig('cell_score.svg', format='svg',bbox_inches='tight')

#Fig S2.D
slice_gene_count = adata.obs.groupby('slice')['n_genes'].apply(lambda x: (x >= 340).mean()).reset_index()
slice_gene_count.columns = ['slice', 'gene_ratio']
plt.figure(figsize=(10, 3))
sns.barplot(data=slice_gene_count, x='slice', y='gene_ratio', color='lightblue')
plt.axhline(y=0.7, color='red', linestyle='--')
plt.title('Proportion of n_genes >= 340 per slice',fontfamily='Arial', fontsize=20)
plt.xlabel('')
plt.ylabel('Proportion',fontfamily='Arial', fontsize=20)
plt.xticks(rotation=45, ha='right',fontsize=12, fontfamily='Arial')
plt.yticks(fontsize=12, fontfamily='Arial')
plt.ylim(0.5, 1.0) 
plt.tight_layout()
plt.savefig('ratio_min_genes_340.svg', format='svg',bbox_inches='tight')

#Fig S2.F
adata=adata[adata.obs['n_genes']>=340]
adata = adata[adata.obs['percent_mt'] <= 0.2]
tmp=pd.read_csv('E8.5_markers.csv')
adata_genes = adata.var_names.tolist()
intersection_genes = list(set(tmp.gene) & set(adata_genes))
min_cells_values = list(range(10, 21))
intersection_ratios = []
for min_cells in min_cells_values:
    adata2=adata.copy()
    sc.pp.filter_genes(adata2,min_cells=min_cells)
    intersection_count = len(adata2.var_names & set(intersection_genes))
    intersection_ratio = intersection_count / 1097
    intersection_ratios.append(intersection_ratio)
plt.plot(min_cells_values, intersection_ratios, marker='o', color='blue')
plt.xlabel('min_cells',fontfamily='Arial', fontsize=20)
plt.ylabel('Proportion',fontfamily='Arial', fontsize=20)
plt.xticks(fontsize=12, fontfamily='Arial')
plt.yticks(fontsize=12, fontfamily='Arial')
plt.title('Ratio of remaining markers vs min_cells',fontfamily='Arial', fontsize=20)
plt.savefig('ratio_min_cells_15.svg', format='svg',bbox_inches='tight')
```

