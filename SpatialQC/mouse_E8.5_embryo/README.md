raw/  
The result of running SpatialQC on the raw data, including the code for obtaining the raw data, the parameters for running SpatialQC, the generated HTML report, and the code for generating Fig S1.  

min_genes_200/  
The result of running SpatialQC on the initially filtered data (filtered mostly non-embryonic tissue cells), including the parameters for running SpatialQC, the generated HTML report, and the code for generating Fig S2.  
The initial filtering uses only the following code:
```python
#python
adata=sc.read_h5ad('E8.5-Atlas_run.h5ad')
sc.pp.filter_cells(adata,min_genes=200)
del adata.obs['n_genes']
adata.write('E8.5-Atlas_run_min200.h5ad')
```
E8.5_markers.csv: marker genes of mouse E8.5 embryo
