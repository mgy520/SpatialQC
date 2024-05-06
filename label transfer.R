###############################
#DLPFC
setwd('/rad/mgy/results/SpatialQC_run_h5ad_dataests/test_h5ad/DLPFC/')
wget https://libd-snrnaseq-pilot.s3.us-east-2.amazonaws.com/SCE_DLPFC-n3_tran-etal.rda
load('SCE_DLPFC-n3_tran-etal.rda')
sc_seu <- as.Seurat(sce.dlpfc.tran, counts = "counts",data = 'logcounts')
sc_seu <- FindVariableFeatures(sc_seu, selection.method = "vst", nfeatures = 3000)

#raw data label transfer
library("anndata")
data <- read_h5ad("dlpfc.h5ad")
data <- CreateSeuratObject(counts = t(data$X), meta.data = data$obs)
data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes)
data <- RunPCA(data)
#anchors <- FindTransferAnchors(reference = sc_seu, query = data)
anchors <- FindTransferAnchors(reference = sc_seu, query = data,k.filter = NA)
predictions.assay <- TransferData(anchorset = anchors, refdata = sc_seu$cellType,
                                weight.reduction = data[["pca"]], dims = 1:30)
data <- AddMetaData(data, metadata = predictions.assay)

#filtered data label transfer
data2 <- read_h5ad("filtered.h5ad")
data2 <- CreateSeuratObject(counts = t(data2$X), meta.data = data2$obs)
data2 <- NormalizeData(data2)
data2 <- FindVariableFeatures(data2, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(data2)
data2 <- ScaleData(data2, features = all.genes)
data2 <- RunPCA(data2)
anchors2 <- FindTransferAnchors(reference = sc_seu, query = data2,k.filter = NA)
predictions.assay2 <- TransferData(anchorset = anchors2, refdata = sc_seu$cellType,
                                  weight.reduction = data2[["pca"]], dims = 1:30)
data2 <- AddMetaData(data2, metadata = predictions.assay2)

values1 <- data$prediction.score.max
values2 <- data2$prediction.score.max
df <- data.frame(value = c(values1, values2), group = factor(c(rep("Raw", length(values1)), rep("Filtered", length(values2))), levels = c("Raw", "Filtered")))




###############################
#MERFISH
wget https://cell2location.cog.sanger.ac.uk/tutorial/mouse_brain_snrna/all_cells_20200625.h5ad
wget https://cell2location.cog.sanger.ac.uk/cellxgene/mouse_brain_snRNA_with_UMAP_cellxgene.h5ad
tmp=sc.read_h5ad('/rad/mgy/results/SpatialQC_run_h5ad_dataests/test_h5ad/MERFISH/mouse_brain_snRNA_with_UMAP_cellxgene.h5ad')
tmp2 = sc.read_h5ad('/rad/mgy/results/SpatialQC_run_h5ad_dataests/test_h5ad/MERFISH/all_cells_20200625.h5ad')
common_cells = list(set(tmp.obs_names) & set(tmp2.obs_names))
tmp2 = tmp2[common_cells, :]
annotation_1_print_a = tmp.obs['annotation_1_print']
tmp2.obs['annotation_1_print'] = tmp2.obs.index.map(annotation_1_print_a)
tmp2.var_names = tmp2.var['SYMBOL']
tmp2.write_h5ad('/rad/mgy/results/SpatialQC_run_h5ad_dataests/test_h5ad/MERFISH/mouse_brain_scrna.h5ad')

setwd('../MERFISH/')
sc_seu <- read_h5ad("mouse_brain_scrna.h5ad")
sc_seu <- CreateSeuratObject(counts = t(sc_seu$X), meta.data = sc_seu$obs,min.cells=3)
sc_seu <- NormalizeData(sc_seu)
sc_seu <- FindVariableFeatures(sc_seu, selection.method = "vst", nfeatures = 3000)

#raw data label transfer
library("anndata")
data <- read_h5ad("vizgen.h5ad")
data <- CreateSeuratObject(counts = t(data$X), meta.data = data$obs)
data <- NormalizeData(data)
all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes)
data <- RunPCA(data,features = all.genes)
anchors <- FindTransferAnchors(reference = sc_seu, query = data)
predictions.assay <- TransferData(anchorset = anchors, refdata = sc_seu$annotation_1_print,
                                   weight.reduction = data[["pca"]], dims = 1:30)
data <- AddMetaData(data, metadata = predictions.assay)

#filtered data label transfer
data2 <- read_h5ad("filtered.h5ad")
data2 <- CreateSeuratObject(counts = t(data2$X), meta.data = data2$obs)
data2 <- NormalizeData(data2)
all.genes <- rownames(data2)
data2 <- ScaleData(data2, features = all.genes)
data2 <- RunPCA(data2,features = all.genes)
anchors2 <- FindTransferAnchors(reference = sc_seu, query = data2)
predictions.assay2 <- TransferData(anchorset = anchors2, refdata = sc_seu$annotation_1_print,
                                    weight.reduction = data2[["pca"]], dims = 1:30)
data2 <- AddMetaData(data2, metadata = predictions.assay2)

values1 <- data$prediction.score.max
values2 <- data2$prediction.score.max
df <- data.frame(value = c(values1, values2), group = factor(c(rep("Raw", length(values1)), rep("Filtered", length(values2))), levels = c("Raw", "Filtered")))





###############################
#xenium
setwd('../xenium/')
wget https://cf.10xgenomics.com/samples/cell-exp/7.0.1/SC3pv3_GEX_Breast_Cancer_DTC_Aggr/SC3pv3_GEX_Breast_Cancer_DTC_Aggr_count_filtered_feature_bc_matrix.tar.gz
tar -zxvf SC3pv3_GEX_Breast_Cancer_DTC_Aggr_count_filtered_feature_bc_matrix.tar.gz
sc_seu <- Read10X('scrna/filtered_feature_bc_matrix/')
sc_seu = CreateSeuratObject(counts = sc_seu)
wget https://cf.10xgenomics.com/samples/cell-exp/7.0.1/SC3pv3_GEX_Breast_Cancer_DTC_Aggr/SC3pv3_GEX_Breast_Cancer_DTC_Aggr_count_analysis.tar.gz
tar -zxvf SC3pv3_GEX_Breast_Cancer_DTC_Aggr_count_analysis.tar.gz
cluster <- read.table('scrna/analysis/clustering/gene_expression_graphclust/clusters.csv',sep = ',',header = TRUE)
sc_seu$cluster <- cluster$Cluster
sc_seu <- NormalizeData(sc_seu)
sc_seu <- FindVariableFeatures(sc_seu, selection.method = "vst", nfeatures = 3000)

#raw data label transfer
library("anndata")
data <- read_h5ad("xenium.h5ad")
data <- CreateSeuratObject(counts = t(data$X), meta.data = data$obs)
data <- NormalizeData(data)
all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes)
data <- RunPCA(data,features = all.genes)
anchors <- FindTransferAnchors(reference = sc_seu, query = data)
sc_seu$cluster <- as.character(sc_seu$cluster)
predictions.assay <- TransferData(anchorset = anchors, refdata = sc_seu$subclass,
                                   weight.reduction = data[["pca"]], dims = 1:30)
data <- AddMetaData(data, metadata = predictions.assay)

#filtered data label transfer
data2 <- read_h5ad("filtered.h5ad")
data2 <- CreateSeuratObject(counts = t(data2$X), meta.data = data2$obs)
data2 <- NormalizeData(data2)
all.genes <- rownames(data2)
data2 <- ScaleData(data2, features = all.genes)
data2 <- RunPCA(data2,features = all.genes)
anchors2 <- FindTransferAnchors(reference = sc_seu, query = data2)
predictions.assay2 <- TransferData(anchorset = anchors2, refdata = sc_seu$cluster,
                                    weight.reduction = data2[["pca"]], dims = 1:30)
data2 <- AddMetaData(data2, metadata = predictions.assay2)

values1 <- data$prediction.score.max
values2 <- data2$prediction.score.max
df <- data.frame(value = c(values1, values2), group = factor(c(rep("Raw", length(values1)), rep("Filtered", length(values2))), levels = c("Raw", "Filtered")))






###############################
#visium_HD
setwd('../visium_HD/')
#scRNA-seq from https://rupress.org/jem/article/220/5/e20221437/213924/Intestinal-cell-type-specific-communication
a <- Read10X(data.dir = 'scrna/')
b <- CreateSeuratObject(counts = a)
b <- subset(b, subset = nFeature_RNA >= 100 & nCount_RNA >= 100)
b[["percent.mt"]] <- PercentageFeatureSet(b, pattern = "^Mt")
b <- subset(b, subset = percent.mt < 2)
b = NormalizeData(b)
b <- FindVariableFeatures(b, selection.method = "vst", nfeatures = 3000)
b <- ScaleData(b, vars.to.regress = 'nFeature_RNA')
b <- RunPCA(b, features = VariableFeatures(object = b))
ElbowPlot(b)
b <- FindNeighbors(b, dims = 1:50)
b <- FindClusters(b,resolution = 0.25)

#raw data label transfer
library("anndata")
data <- read_h5ad("dlpfc.hd5ad")
data <- CreateSeuratObject(counts = t(data$X), meta.data = data$obs)
data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes)
data <- RunPCA(data)
anchors <- FindTransferAnchors(reference = b, query = data)
predictions.assay <- TransferData(anchorset = anchors, refdata = b$seurat_clusters,
                                   weight.reduction = data[["pca"]], dims = 1:30)
data <- AddMetaData(data, metadata = predictions.assay)

#filtered data label transfer
data2 <- read_h5ad("filtered.h5ad")
data2 <- CreateSeuratObject(counts = t(data2$X), meta.data = data2$obs)
data2 <- NormalizeData(data2)
data2 <- FindVariableFeatures(data2, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(data2)
data2 <- ScaleData(data2, features = all.genes)
data2 <- RunPCA(data2)
anchors2 <- FindTransferAnchors(reference = b, query = data2)
predictions.assay2 <- TransferData(anchorset = anchors2, refdata = b$seurat_clusters,
                                   weight.reduction = data2[["pca"]], dims = 1:30)
data2 <- AddMetaData(data2, metadata = predictions.assay2)

values1 <- data$prediction.score.max
values2 <- data2$prediction.score.max
df <- data.frame(value = c(values1, values2), group = factor(c(rep("Raw", length(values1)), rep("Filtered", length(values2))), levels = c("Raw", "Filtered")))





###############################
#Stereo-seq
setwd('../Drosophila/')
b <- readRDS('scrna1/14_16_finished_processing.Rds')
b <- SCTransform(b)

#raw data label transfer
library("anndata")
data <- read_h5ad("Drosophila_embryo_E14-16h_spatialqc_test.h5ad")
data <- CreateSeuratObject(counts = t(data$X), meta.data = data$obs)
data <- SCTransform(data,min_cells=1)
data <- RunPCA(data,features = VariableFeatures(object = data))
sc_seu = b
anchors <- FindTransferAnchors(reference = sc_seu, query = data,normalization.method = 'SCT')
predictions.assay <- TransferData(anchorset = anchors, refdata = sc_seu$seurat_clusters,
                                   weight.reduction = data[["pca"]], dims = 1:30)
data <- AddMetaData(data, metadata = predictions.assay)

#filtered data label transfer
data2 <- read_h5ad("filtered.h5ad")
data2 <- CreateSeuratObject(counts = t(data2$X), meta.data = data2$obs)
data2 <- SCTransform(data2,min_cells=1)
data2 <- RunPCA(data2)
anchors2 <- FindTransferAnchors(reference = sc_seu, query = data2,normalization.method = 'SCT')
predictions.assay2 <- TransferData(anchorset = anchors2, refdata = sc_seu$seurat_clusters,
                                    weight.reduction = data2[["pca"]], dims = 1:30)
data2 <- AddMetaData(data2, metadata = predictions.assay2)

values1 <- data$prediction.score.max
values2 <- data2$prediction.score.max
df <- data.frame(value = c(values1, values2), group = factor(c(rep("Raw", length(values1)), rep("Filtered", length(values2))), levels = c("Raw", "Filtered")))






###############################
#Slide-seq
setwd('/rad/mgy/results/NG_2023_spatialqc')
a <- readRDS('/rad/share/stereo-seq/embryo_rds/mouse_gastrulation_recluster.Rds')
a2 <- a[,a$stage=='E8.5']
a2 <- a2[,a2$stripped==FALSE]
a2 <- a2[,a2$doublet==FALSE]
rm(a)

#raw data label transfer
library("anndata")
data <- read_h5ad("E8.5-Atlas_run_min200.h5ad")
data <- CreateSeuratObject(counts = t(data$X), meta.data = data$obs)
data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
data <- ScaleData(data, features = VariableFeatures(object = data))
data <- RunPCA(data,features = VariableFeatures(object = data))
sc_seu = a2
anchors <- FindTransferAnchors(reference = sc_seu, query = data)
predictions.assay <- TransferData(anchorset = anchors, refdata = sc_seu$celltype,
                                   weight.reduction = data[["pca"]], dims = 1:30)
data <- AddMetaData(data, metadata = predictions.assay)

#filtered data label transfer
data2 <- read_h5ad("filtered.h5ad")
data2 <- CreateSeuratObject(counts = t(data2$X), meta.data = data2$obs)
data2 <- NormalizeData(data2)
data2 <- FindVariableFeatures(data2, selection.method = "vst", nfeatures = 2000)
data2 <- ScaleData(data2, features = VariableFeatures(object = data2))
data2 <- RunPCA(data2,features = VariableFeatures(object = data2))
anchors2 <- FindTransferAnchors(reference = sc_seu, query = data2)
predictions.assay2 <- TransferData(anchorset = anchors2, refdata = sc_seu$celltype,
                                    weight.reduction = data2[["pca"]], dims = 1:30)
data2 <- AddMetaData(data2, metadata = predictions.assay2)

values1 <- data$prediction.score.max
values2 <- data2$prediction.score.max
df <- data.frame(value = c(values1, values2), group = factor(c(rep("Raw", length(values1)), rep("Filtered", length(values2))), levels = c("Raw", "Filtered")))