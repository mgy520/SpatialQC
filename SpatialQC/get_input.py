import anndata2ri
from rpy2.robjects import r

anndata2ri.activate()


def convert_to_anndata(input1):
    r.assign("input", input1)
    r(
    '''
    library(SpatialExperiment)
    library(Seurat)
    library(SingleCellExperiment)
    input2 = readRDS(input)
    input_class = class(input2)
    if (input_class == "Seurat") {
        if (length(input2@images) > 0) {
            coords = GetTissueCoordinates(input2)
            colnames(coords)[1:2] <- c("row", "col")
            input2@meta.data = cbind(input2@meta.data, coords)
        }
        input2 = as.SingleCellExperiment(input2)
    } else if (input_class == "SpatialExperiment") {
        colnames(spatialCoords(input2)) <- c('row','col')
        colData(input2) <- cbind(colData(input2), spatialCoords(input2))
    }      
    adata = as(input2, "SingleCellExperiment")
    '''
    )
    adata = r['adata']
    adata.obsm['spatial'] = adata.obs[['row', 'col']]

    return adata
