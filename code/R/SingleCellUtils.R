library(Seurat)
library(patchwork)
library(enrichplot)
library(org.Hs.eg.db)
library(DOSE)

FilterDegs <- function(markers, seurat_object, n_genes = 6){
  # Initialize an empty list for storing DEGs (differentially expressed genes)
  deg_list <- list()
  cluster_num <- length(levels(seurat_object)) -1 
  # Iterate through clusters (0 to 4)
  for (i in 0:cluster_num) {
    # Filter markers based on avg_log2FC and cluster
    degs <- markers %>%
      group_by(cluster) %>%
      filter(avg_log2FC > 1, cluster == i) %>%
      arrange(p_val_adj) %>%
      slice_head(n = n_genes) %>%
      pull(gene)
    # Check if any DEGs were found
    degs = sort(degs, decreasing = TRUE)
    deg_list[[i+1]] <- degs
  }
  return(deg_list)
}


ClusterComp <- function (seurat_object){
  # Get the cluster identities and original identities
  Idents(seurat_object) <- 'RNA_snn_res.0.4'
  cluster_ids <- Idents(seurat_object)
  orig_ids <- seurat_object$orig.ident
  
  # Create a data frame
  df <- data.frame(Cluster = cluster_ids, OrigIdent = orig_ids)
  
  # Count the number of cells from each orig.ident in each cluster
  df_counts <- table(df)
  
  # Convert the table to a data frame
  df_counts <- as.data.frame(df_counts)
  
  # Rename the columns
  colnames(df_counts) <- c("Cluster", "OrigIdent", "Count")
  
  # Create the fraction column
  df_counts <- df_counts %>%
    group_by(Cluster) %>%
    mutate(TotalCount = sum(Count), Fraction = Count / TotalCount) %>%
    ungroup() # Ungroup to avoid issues in further data manipulations
  
  count_plot <- ggplot(df_counts, aes(x = Cluster, y = Count, fill = OrigIdent)) +
    geom_bar(stat = "identity") +
    theme(legend.position = "none") +
    labs(x = "Cluster", y = "Cluster Composition (Cell Count)", fill = "Sample", title = "Sample Composition of clusters")
  
  fraction_plot <- ggplot(df_counts, aes(x = Cluster, y = Fraction, fill = OrigIdent)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(x = "Cluster", y = "Cluster Composition (%)", fill = "Sample", title = "Sample Composition of clusters")
  
  count_plot + fraction_plot
}

GO_KEGG <- function(seurat_object){
  

  deg.ls <- split(rownames(seurat_object), f = seurat_object$RNA_snn_res.0.4)
  
  geneid.ls <- deg.ls %>% map(~{
    
    # here for macaque
    gene.df <- select(org.Hs.eg.db,
                      keys = .x,
                      columns = c("ENTREZID", "SYMBOL"),
                      keytype = "SYMBOL")
    
    gene <- gene.df$ENTREZID
    gene <- gene[which(!is.na(gene))]
    gene <- unique(gene)
    
  })
  
  gene.ls <- geneid.ls
  
  # her mcc for macaque
  compKEGG <- compareCluster(geneCluster   = gene.ls,
                             fun           = "enrichKEGG",
                             pvalueCutoff  = 0.05,
                             pAdjustMethod = "BH", 
                             organism = "hsa")
  
  compGO <- compareCluster(geneCluster   = gene.ls,
                           fun           = "enrichGO",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH", 
                           OrgDb = org.Hs.eg.db, 
                           ont = 'BP')
  
  #  compPathway <- compareCluster(geneCluster   = gene.ls,
  #                                fun           = "enrichPathway",
  #                                pvalueCutoff  = 0.05,
  #                                pAdjustMethod = "BH")
  
  ## dot plot
  g1 <- dotplot(compGO, showCategory = 10, title = "GO Enrichment Analysis")
  #g2 <- dotplot(compPathway, showCategory = 10, title = "REACTOME Pathway Enrichment Analysis")
  g3 <- dotplot(compKEGG, showCategory = 10, title = "KEGG Pathway Enrichment Analysis")
  plots <- c()
  plots[[1]] <- g1
  plots[[2]] <- g3
  return (plots)
}

CreateSeuratList <- function(file_paths){
  # Initialize an empty list to store the Seurat objects for each sample
  seurat_list <- list()
  
  for (i in seq_along(file_paths)){
    seurat_data <- Read10X(data.dir = file_paths[[i]])
    sample_name <- strsplit(file_paths[[i]],"/")[[1]][[4]]
    seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                     min.features = 750, 
                                     min.cells = 3,
                                     project = sample_name)
    seurat_list[[i]] <- seurat_obj
  }
  
  merged_seurat <- merge(x = seurat_list[[1]], 
                         y = seurat_list[-1])
  
  # Remove objects to clear up memory
  rm(seurat_data)
  rm(seurat_list)
  rm(seurat_obj)
  
  # Compute percent mito ratio
  merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
  merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100
  
  
  # Add number of genes per UMI for each cell to metadata
  merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)
  
  return (merged_seurat)
}