## **Main arguments**

`--input`  None  
The path to your .h5ad or .rds file. Supports object: anndata, Seurat, SingleCellExperiment, SpatialExperiment.
Note that the spatial coordinates of anndata objects are stored in obsm['spatial']. If your Seurat object does not have an images slot, the row and col coordinates need to be provided in meta.data. SingleCellExperiment must provide row and col coordinates.


`--platform`  None  
Platform of your spatial transcriptomics data. Visium, MERFISH, Stereo-seq, etc.  
Visium: Same as `--doublet False --n 0.95 --l 0.99 --s 4 --min_genes_list 0 200 400 600 800 1000 1200
--min_genes_list2 0 200 400 600 800 1000 1200 --min_cells_list 1 2 3 4 5`. Or set the parameter to 'Slide-seq', 'ST', 'DBiT-seq'.  
MERFISH: Same as `--doublet True --n 0.9 --min_cells 1 --s 3 --min_genes_list 0 10 20 30 40 50 
--min_genes_list2 0 10 20 30 40 50 --min_cells_list 1 --s2 0 0 1`. Or set the parameter to 'Xenium', 'CosMx', 'HybISS'.  
Stereo-seq: Same as `--doublet True --n 0.7 --l 0.99 --s4 --min_genes_list 0 200 400 600 800 1000 1200
--min_genes_list2 0 200 400 600 800 1000 1200 --min_cells_list 0 10 20 30 40`. Or set the parameter to 'Seq-scope', 'Pixel-seq', 'HDST', 'Visium HD'.

`--slice_number`  multiple  
The number of slices of .h5ad provided. multiple or 1. 

`--slice`   id  
The name that represents the slice identifier in `anndata.obs`. If there is only one slice, ignore this parameter.

`--mito`  'Mt-'  
The pattern of mitochondrial genes. 

`--ribo`  'Rps, Rpl'  
The pattern of ribosome genes. 

`--hemo`  'Hbb, Hba'  
The pattern of hemoglobin genes. 



## **Marker genes arguments**

`--markers`  None  
The path to your marker genes .csv file. 
You can obtain genes based on prior knowledge or DEGs 
from scRNA-seq of the same tissue.  
We also provided mouse and human marker genes obtained from the 
cellmarker2.0 database. If you don't have a suitable markers file, you can specify the parameters `--species`, `--tissue_class`,
`--tissue_type`, `--cancer_type`.  
If markers are not provided, set to False.

`--species`  None  
The species of your sample. Human or Mouse. 

`--tissue_class`  None  
The tissue class of your sample. View options in [CellMarker2.0](http://117.50.127.228/CellMarker/index.html). 

`--tissue_type`  None  
The tissue type of your sample. View options in [CellMarker2.0](http://117.50.127.228/CellMarker/index.html). 

`--cancer_type`  Normal  
The cancer type of your sample. View options in [CellMarker2.0](http://117.50.127.228/CellMarker/index.html). 



## **Filter arguments**

`--f`  True  
Whether to filter .h5ad file. if False, generate only HTML report.

`--s`  5  
Sections with a median score less than s will be removed.

`--n`  0.7  
Determine the value of min_genes to ensure that the valid cell ratio is greater
than `--n`. min_genes is adjusted to the nearest multiple of 10. If min_genes
is already divisible by 10, it remains unchanged. Otherwise, min_genes is 
rounded down to the nearest multiple of 10.

`--min_genes`  None  
Provide your min_genes, otherwise determined by `--n`.

`--l`  0.99  

After filtering cells, determine the value of min_cells to ensure that the proportion of marker genes is greater than `--l` among the remaining detected markers.

`--min_cells`  None  
Provide your min_cells, otherwise determined by `--l`.

`--mito_percent`  0.1  
Filter cells with mitochondrial proportion higher than `--mito_percent`.



## **Cell score arguments**

`--s1`  -1 0  
percent.mt_score of the cell.   
percent_mito > `--mito_percent`: -1  
percent_mito <= `--mito_percent`: 0

 `--s2`  0.8 0 1  
log10GenesPerUMI_score of the cell.  
log10GenesPerUMI < 0.8: 0  
log10GenesPerUMI >= 0.8: 1

`--s3`  0.2 0.5 0 1 2  
n_genes_score of the cell.  
n_genes ranking below 20th percentile: 0  
between 20th and 50th percentile: 1  
above 50th percentile: 2

`--s4`  0.2 0.5 0.8 0 1 2 3  
markerDetectionRatio_score of the cell.  
markerDetectionRatio ranking below 20th percentile: 0  
between 20th and 50th percentile: 1  
between 50th and 80th percentile: 2  
above 80th percentile: 3

`--s5`  0.2 0.5 0.8 0 1 2 3]  
markerProportion_score of the cell.  
markerProportion_score ranking below 20th percentile: 0  
 between 20th and 50th percentile: 1  
between 50th and 80th percentile: 2  
above 80th percentile: 3

`--s6`  0.2 0.5 0.8 0 1 2 3  
markerCountsRatio_score of the cell.  
markerCountsRatio_score ranking below 20th percentile: 0  
between 20th and 50th percentile: 1  
between 50th and 80th percentile: 2  
above 80th percentile: 3

`--s7`  -4 0  
doublet_score of the cell.  
doublet cells: -4  
not doublet cells: 0

`--s8`  0.2 0.5 0 1 2  
n_counts_score.  
n_counts ranking below 20th percentile: 0  
between 20th and 50th percentile: 1  
above 50th percentile: 2



## **Options**

`--bin_value`  100  
Values of n_genes bin intervals applied to the HTML button: Marker Proportion.

`--min_genes_list`  0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100  
used for HTML buttons: Cell Number Post Filter, Markers Proportion Post Filter, Markers Detected Post Filter. 

`--min_genes_list2`  0, 100, 200, 300, 400, 500, 600, 700  
used for HTML button: Valid Cell Post min_genes.

`--min_cells_list`  3, 10, 20  
used for html buttons: Markers Proportion Post Filter, Markers Detected Post Filter.

`--output`  ./  
output directory.

`--o1`  report.html  
The filename for the output of html report.

`--o2`  filtered.h5ad  
The filename for the output of filtered .h5ad.

`--j`  8  
The maximum number of concurrently running jobs. If set to 1, parallelism is not used. If set to -1, all CPUs are used. For n_jobs less than -1, (n_cpus + 1 + n_jobs) CPUs are used.