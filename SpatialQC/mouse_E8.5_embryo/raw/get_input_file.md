### 1. Obtain the input mouse E8.5 embryo 'anndata' format data

The raw data was downloaded from GSE197353

```R
setwd('NG_2023')
E85_pucks = c("201104_07","201104_08","201104_09","201104_12","201104_13",
              "201104_14","201104_16","201104_17","201104_18","201104_19",
              "201104_22","201104_23","201104_24","201104_26","201104_27",
              "201104_28","201104_29")
srt_list = list()
for (sel_puck in E85_pucks) {
  print(sel_puck)
  pos = read.csv(paste0("GSE197353_RAW/",sel_puck,"_matched_bead_locations.txt.gz"),sep="\t",header = F)
  cnt = read.csv(paste0("GSE197353_RAW/",sel_puck,".digital_expression.txt.gz"),sep="\t",row.names = 1)
  
  pos = pos[,2:3]
  rownames(pos) = colnames(cnt)
  colnames(pos) = c("xcoord","ycoord")
  srt = CreateSeuratObject(counts = cnt,meta.data = pos,min.features=0,project = sel_puck)
  srt$stage = "E8.5"
  srt_list[[sel_puck]] = srt
}

#merge
srt = merge(srt_list[[1]],y=srt_list[-1],add.cell.ids=names(srt_list))
emb = srt@meta.data[,c("xcoord","ycoord")]
colnames(emb) = c("s_1","s_2")
srt[["spatial"]] <- CreateDimReducObject(embeddings = as.matrix(emb), key = "s_")
MuDataSeurat::WriteH5AD(srt, "E8.5-Atlas_raw.h5ad")
```

Check that python works
```python
import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read('E8.5-Atlas_raw.h5ad')
#add 'slice' key to adata.obs
def extract_number(string):
    return string.split('_')[-1]

def add_slice(orig_ident):
    number = extract_number(orig_ident)
    return f"sagittal_{number}"

adata.obs['slice'] = adata.obs['orig.ident'].apply(add_slice)
xcoord = adata.obs['xcoord']
ycoord = adata.obs['ycoord']
adata.obsm['spatial'] = np.column_stack([xcoord, ycoord])
del adata.obs['nCount_RNA']
del adata.obs['nFeature_RNA']
del adata.layers
adata.write('E8.5-Atlas_run.h5ad')
```

### 2. Get marker genes

```R
a <- readRDS('mouse_gastrulation_recluster.Rds')
a2 <- a[,a$stage=='E8.5']
a2 <- a2[,a2$stripped==FALSE]
a2 <- a2[,a2$doublet==FALSE]
markers <- FindAllMarkers(a2, logfc.threshold = 0.5, only.pos = TRUE)
top100 <- markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
top100 <- top100[,'gene']
write.table(top100,file = 'E8.5_markers.csv',quote = F,row.names = F)
```

