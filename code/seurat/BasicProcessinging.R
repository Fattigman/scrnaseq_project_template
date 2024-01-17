# Load libraries
library(Seurat)
library(dplyr)
source("code/seurat/ClusterComp.R")
# List of file paths to the Cell Ranger output directories
# Get all directories
all_dirs <- list.dirs(path = "workspace/data/raw_internal", full.names = TRUE)
file_paths <- paste(all_dirs, "/outs/filtered_feature_bc_matrix", sep = "")

# Initialize an empty list to store the Seurat objects for each sample
seurat_list <- list()

for (i in seq_along(file_paths)){
  seurat_data <- Read10X(data.dir = file_paths[[i]])
  sample_name <- strsplit(file_paths[[i]],"/")[[1]][[4]]
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 500, 
                                   min.cells = 3,
                                   project = sample_name)
  seurat_list[[i]] <- seurat_obj
}

merged_seurat <- merge(x = seurat_list[[1]],
                       y = seurat_list[-1],
                       add.cell.id = c("jp_baseline",
                                       "JP_w4", "kp_baseline", "kp_w4"))

# Remove objects to clear up memory
rm(seurat_data)
rm(seurat_list)
rm(seurat_obj)

# Uncomment and modify to add metadata
#merged_seurat <- merged_seurat %>%
#  AddMetaData(metadata = ifelse(grepl("JP", merged_seurat@meta.data$orig.ident), "JP", "KP"), col.name = "patient")

# Add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

# Compute percent mito ratio
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100

# Filter out low quality cells using selected thresholds - these will change with experiment
filtered_seurat <- subset(x = merged_seurat, 
                          subset= (
                            nCount_RNA >= 500) & 
                            (nFeature_RNA >= 250) & 
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.20))


# Remove objects to clear up memory
rm(merged_seurat)

# Normalize the data
seurat_phase <- NormalizeData(filtered_seurat)
seurat_phase <- JoinLayers(seurat_phase)
# Score cells for cell cycle
seurat_phase <- CellCycleScoring(seurat_phase,
                                 g2m.features = cc.genes$g2m.genes,
                                 s.features = cc.genes$s.genes,
                                 set.ident = TRUE)

# Identify the most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase,
                                     selection.method = "vst",
                                     nfeatures = 2000,
                                     verbose = FALSE)

# Scale the counts
seurat_phase <- ScaleData(seurat_phase)

# Perform PCA
seurat_phase <- RunPCA(seurat_phase)

# Turn mitoRatio into categorical factor vector based on quartile values
seurat_phase@meta.data$mitoFr <- cut(seurat_phase@meta.data$mitoRatio, 
                                     breaks=c(-Inf, 0.0144, 0.0199, 0.0267, Inf),
                                     labels=c("Low","Medium","Medium high", "High"))

# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
split_seurat <- SplitObject(seurat_phase, split.by = "orig.ident")

options(future.globals.maxSize = 4000 * 1024^2)

for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vst.flavor = "v2")
}

# Run UMAP
seurat_phase <- RunUMAP(seurat_phase,
                        dims = 1:15,reduction = "pca")
# Plot UMAP
DimPlot(seurat_phase, group.by = "patient")

# Determine the K-nearest neighbor graph
seurat_phase <- FindNeighbors(object = seurat_phase,
                              dims = 1:10)

seurat_phase <- FindClusters(object = seurat_phase,
                             resolution = c(0.4))

Idents(seurat_phase) <- seurat_phase@meta.data$RNA_snn_res.0.4

DimPlot(seurat_phase,reduction = "umap", label = TRUE)

# Check the cluster composition
ClusterComp(seurat_phase)


DimPlot(seurat_phase, group.by = c("patient", "RNA_snn_res.0.4"), label = TRUE)

# Automatic cell type annotater
library(openxlsx)
library(HGNChelper)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# get cell-type-specific gene sets from our in-built database (DB)
gs_list = gene_sets_prepare("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx", "Immune system") # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

es.max = sctype_score(scRNAseqData = seurat_phase[['RNA']]$scale.data, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(seurat_phase@meta.data$RNA_snn_res.0.4), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(seurat_phase@meta.data[seurat_phase@meta.data$RNA_snn_res.0.4==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_phase@meta.data$RNA_snn_res.0.4==cl)), 10)
}))

sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

# Plot the clusters on umap
seurat_phase@meta.data$celltype = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  seurat_phase@meta.data$celltype[seurat_phase@meta.data$RNA_snn_res.0.4 == j] = as.character(cl_type$type[1])
}

DimPlot(seurat_phase,
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  group.by = c("celltype", "patient")
) + ggtitle("Cell type classification")