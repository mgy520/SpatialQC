### The code for generating Fig S1

```python
#import the required packages 
import scanpy as sc
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#load raw data
adata = sc.read('E8.5-Atlas_run.h5ad')
sc.pp.filter_cells(adata,min_genes=0)

#Fig S1.A
plt.figure()
tmp=adata.obs['n_genes']
sns.kdeplot(tmp, fill=False,log_scale=True,color='blue')
plt.axvline(x=200, color='red', linestyle='-')
plt.xlabel('n_genes', fontfamily='Arial', fontsize=20)
plt.ylabel('log10 cell density', fontfamily='Arial', fontsize=20)
plt.xticks(fontsize=12, fontfamily='Arial')
plt.yticks(fontsize=12, fontfamily='Arial')
plt.text(250, 0.55, 'n_genes=200', color='red', fontsize=10, fontfamily='Arial')
plt.savefig('n_genes_density_plot.svg', format='svg',bbox_inches='tight')

#Fig S1.B
sc.pp.filter_cells(adata,min_counts=0)
plt.figure()
tmp=adata.obs['n_counts']
sns.kdeplot(tmp, fill=False,log_scale=True,color='blue')
plt.axvline(x=500, color='red', linestyle='-')
plt.xlabel('n_counts', fontfamily='Arial', fontsize=20)
plt.ylabel('log10 cell density', fontfamily='Arial', fontsize=20)
plt.xticks(fontsize=12, fontfamily='Arial')
plt.yticks(fontsize=12, fontfamily='Arial')
plt.text(550, 0.55, 'n_counts=500', color='red', fontsize=10, fontfamily='Arial')
plt.savefig('n_counts_density_plot.svg', format='svg',bbox_inches='tight')

#Fig S1.C
plt.figure()
tmp = adata.obs['n_genes']
sns.violinplot(y=tmp,facecolor=(220/255, 65/255, 80/255, 0.5))
plt.title('Violin Plot of n_genes',fontfamily='Arial', fontsize=20)
plt.ylabel('n_genes', fontfamily='Arial', fontsize=20)
def format_func(value, tick_number):
    if value >= 1000:
        value = int(value / 1000)
        return f'{value}k'
    else:
        return value

plt.gca().yaxis.set_major_formatter(FuncFormatter(format_func))
plt.xticks(fontsize=12, fontfamily='Arial')
plt.yticks(fontsize=12, fontfamily='Arial')
plt.savefig('n_genes_box_plot.svg', format='svg',bbox_inches='tight')
plt.show()

#Fig S1.D
plt.figure()
tmp = adata.obs['n_counts']
sns.violinplot(y=tmp,facecolor=(220/255, 65/255, 80/255, 0.5))
plt.title('Violin Plot of n_counts',fontfamily='Arial', fontsize=20)
plt.ylabel('n_counts', fontfamily='Arial', fontsize=20)
def format_func(value, tick_number):
    if value >= 1000:
        value = int(value / 1000)
        return f'{value}k'
    else:
        return value

plt.gca().yaxis.set_major_formatter(FuncFormatter(format_func))
plt.xticks(fontsize=12, fontfamily='Arial')
plt.yticks(fontsize=12, fontfamily='Arial')
plt.savefig('n_counts_box_plot.svg', format='svg',bbox_inches='tight')

#Fig S1.E
tmp = adata.obs[adata.obs['slice'] == 'sagittal_07']
xcoord = tmp['xcoord']
ycoord = tmp['ycoord']
n_genes = tmp['n_genes']
plt.figure(figsize=(5, 3))
plt.scatter(xcoord[n_genes < 200], ycoord[n_genes < 200], color='red', label='n_genes<200',s=0.5)
plt.scatter(xcoord[n_genes >= 200], ycoord[n_genes >= 200], color='green', label='n_genes>=200 ',s=0.5)
plt.axis('equal')
plt.axis('off')
plt.savefig('s07_ngenes_200divide.svg', format='svg',bbox_inches='tight')
plt.savefig('s07_ngenes_200divide.tif', format='tif', dpi=350)

#Fig S1.F
tmp = adata.obs[adata.obs['slice'] == 'sagittal_28']
xcoord = tmp['xcoord']
ycoord = tmp['ycoord']
n_genes = tmp['n_genes']
plt.figure(figsize=(5, 3))
plt.scatter(xcoord[n_genes < 200], ycoord[n_genes < 200], color='red', label='n_genes<200',s=0.5)
plt.scatter(xcoord[n_genes >= 200], ycoord[n_genes >= 200], color='green', label='n_genes>=200 ',s=0.5)
plt.axis('equal')
plt.axis('off')
plt.savefig('s28_ngenes_200divide.svg', format='svg',bbox_inches='tight')
plt.savefig('s28_ngenes_200divide.tif', format='tif', dpi=350)
```

