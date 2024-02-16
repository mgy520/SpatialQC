SpatialQC is a spatial transcriptomic data quality control tool. Available for Linux, Windows and macOS. Run SpatialQC and you'll get an HTML interactive report and clean data.
## Install
```
pip install SpatialQC
```
## Help Documents
```
SpatialQC -h
```
See [https://mgy520.github.io/SpatialQC/](https://mgy520.github.io/SpatialQC/) for detailed help.
## Run command:
```
SpatialQC --adata your.h5ad --markers your_genes.csv
```
## You will get a clean .h5ad file and an HTML interactive report
![image](https://github.com/mgy520/SpatialQC/blob/main/SpatialQC/paper_example/report.png)
Sample reports are available in the SpatialQC/paper_example directory.
