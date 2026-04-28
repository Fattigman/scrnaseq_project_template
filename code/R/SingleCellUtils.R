library(purrr)
library(Seurat)
library(patchwork)
library(sctransform)
library(AUCell)
library(ggplot2)
library(glue)
library(dplyr)
library(ggnewscale)
library(tidyverse)
#library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyr)
library(enrichR)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)

options(omnipath.logfile = "none")


PrepMarkers <- function(markers){
  markers <- markers %>% mutate(pct.diff = pct.1-pct.2) %>% arrange(desc(pct.diff))
  return (markers)
}
FeatureLoadings <- function(seurat_object, gene_list){
  feature.loadings <- seurat_object@reductions$pca@feature.loadings
  genes <- gene_list[gene_list %in% rownames(feature.loadings)]
  if (length(genes) == 0){
    print("No matching genes in feature loading")
  }else if(length(genes) == 1){
    pca_loadings_sums <- seurat_object@reductions$pca@feature.loadings[genes, ]
    pca_df <- data.frame(PC = names(pca_loadings_sums), LoadingSum = pca_loadings_sums)
    pca_df <- pca_df %>%
      mutate(PC_num = as.numeric(gsub("PC_", "", PC))) %>%  # Extract number from "PC_X"
      arrange(PC_num)  # Sort by numeric PC number
    pca_df$PC <- factor(pca_df$PC, levels = pca_df$PC)
    ggplot(pca_df, aes(x = PC, y = LoadingSum)) +
      geom_bar(stat = "identity") +
      theme_minimal() +
      labs(title = "Sum of PCA Loadings for Selected Genes",
           x = "Principal Component",
           y = "Loading Sum") + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
  }else{
    pca_loadings_sums <- colSums(seurat_object@reductions$pca@feature.loadings[genes, ])
    pca_df <- data.frame(PC = names(pca_loadings_sums), LoadingSum = pca_loadings_sums)
    pca_df <- pca_df %>%
      mutate(PC_num = as.numeric(gsub("PC_", "", PC))) %>%  # Extract number from "PC_X"
      arrange(PC_num)  # Sort by numeric PC number
    pca_df$PC <- factor(pca_df$PC, levels = pca_df$PC)
    ggplot(pca_df, aes(x = PC, y = LoadingSum)) +
      geom_bar(stat = "identity") +
      theme_minimal() +
      labs(title = "Sum of PCA Loadings for Selected Genes",
           x = "Principal Component",
           y = "Loading Sum") + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
  }
}
FilterDegs <- function(markers, seurat_object, n_genes = 6){
  # Initialize an empty list for storing DEGs (differentially expressed genes)
  deg_list <- list()
  if ("cluster" %in% names(markers)){
    clusters <- levels(seurat_object)
    # Iterate through clusters (0 to n)
    flag <- 1
    for (i in clusters) {
      # Filter markers based on avg_log2FC and cluster
      degs <- markers %>%
        group_by(cluster) %>%
        filter(abs(avg_log2FC) > 1, cluster == i) %>%
        arrange(p_val_adj) 
      # Check if any DEGs were found
      #degs = sort(degs, decreasing = TRUE)
      deg_list[[flag]] <- head(degs[degs$avg_log2FC > 1,]$gene, n_genes)
      flag <- flag + 1
    }
  }else {
    
    degs <- markers %>%
      filter(abs(avg_log2FC) > 1) %>%
      arrange(p_val_adj) %>%
      slice_head(n = n_genes) %>%
      arrange(avg_log2FC)
    
    # Check if any DEGs were found
    deg_list[[1]] <- rownames(degs[degs$avg_log2FC > 1,])
    deg_list[[2]] <- rownames(degs[degs$avg_log2FC < -1,])
  }
  
  return(deg_list)
}


ClusterComp <-  function(seurat_object, meta_data= "orig.ident", idents = NULL, cluster_name = "seurat_clusters", reduced = FALSE){
  if (reduced) {
    least_cells_sample <- rownames(as.matrix(table(seurat_object$orig.ident)[table(seurat_object$orig.ident) == min(table(seurat_object$orig.ident))]))
    least_cells <- table(seurat_object$orig.ident)[table(seurat_object$orig.ident) == min(table(seurat_object$orig.ident))][1]
    cell_list <- c()
    for (sample in unique(seurat_object$orig.ident)) {
      # Create a subset of the data for the current sample
      if (sample == least_cells_sample){
        next
      }
      subset_sample <- seurat_object[,seurat_object$orig.ident == sample]
      
      # Sample 6413 column names from the subset
      sampled_columns <- sample(colnames(subset_sample), least_cells)
      
      # Append the sampled columns to the cell_list
      cell_list <- append(cell_list, sampled_columns)
    }
    gc()
    seurat_object <- seurat_object[,cell_list]
  }
  # Get the cluster identities and original identities
  cluster_ids <- seurat_object[[cluster_name]]
  orig_ids <- seurat_object@meta.data[[meta_data]]
  
  # Create a data frame
  df <- data.frame(Cluster = cluster_ids, OrigIdent = orig_ids)
  
  # Count the number of cells from each orig.ident in each cluster
  df_counts <- table(df)
  
  # Convert the table to a data frame
  df_counts <- as.data.frame(df_counts)
  
  # Rename the columns
  colnames(df_counts) <- c("Cluster", "OrigIdent", "Count")
  
  # Filter out cells if idents is specified
  if (!is.null(idents)) {
    df_counts <- df_counts %>%
      filter(Cluster %in% idents)
  }
  
  # Create the fraction column
  df_counts <- df_counts %>%
    group_by(Cluster) %>%
    mutate(TotalCount = sum(Count), Fraction = Count / TotalCount) %>%
    ungroup() # Ungroup to avoid issues in further data manipulations
  
  # Create the count plot
  count_plot <- ggplot(df_counts, aes(x = Cluster, y = Count, fill = OrigIdent)) +
    geom_bar(stat = "identity") +
    theme(axis.text = element_text(angle = 90)) +
    labs(x = "Cluster", y = "Cluster composition (Cell Count)", fill = "Sample", angle = 90,title = glue("{meta_data} Cluster composition (Cell Count)"))
  
  # Annotate the bars with counts
  count_plot <- count_plot
  
  # Create the fraction plot
  fraction_plot <- ggplot(df_counts, aes(x = Cluster, y = Fraction, fill = OrigIdent)) +
    geom_bar(stat = "identity") +
    theme(axis.text = element_text(angle = 90)) +
    labs(x = "Cluster", y = "Cluster Composition (%)", fill = "Sample", title = glue("{meta_data} cluster composition (%)"))
  
  count_plot + fraction_plot
}

ClusterComp2 <- function(seurat_object, 
                         x = "Group", 
                         hue = "patient_id", 
                         split = "seurat_clusters",
                         scales = "free",
                         reduce = TRUE,
                         fraction = TRUE){
  meta_df <- seurat_object@meta.data
  meta_df['X'] <- meta_df[[x]]
  meta_df['Hue'] <- meta_df[[hue]]
  meta_df['Split'] <- meta_df[[split]]
  
  if (reduce){
    # Get total counts per X
    total_counts <- meta_df %>% dplyr::count(X, name = "n")
    
    # Get the minimum count across X
    min_n <- min(total_counts$n)
    
    # Calculate scaling factor per X
    scaling_factors <- total_counts %>%
      mutate(scaling_factor = min_n / n) %>%
      select(X, scaling_factor)
    
    # Get counts per X, Split, Hue
    counts_df <- meta_df %>%
      dplyr::count(X, Split, Hue, name = "count")
    
    # Adjust counts and compute fractions within each X
    data <- counts_df %>%
      left_join(scaling_factors, by = "X") %>%
      mutate(scaled_count = count * scaling_factor) %>%
      group_by(X) %>%
      mutate(Fraction = scaled_count / sum(scaled_count)) %>%
      ungroup()
  } else {
    freq_df <- meta_df %>% dplyr::count(X) %>%
      dplyr::summarize(sum_of_cells = sum(n))
    
    data <- meta_df %>%
      dplyr::count(X, Split, Hue, name = "count") %>%
      tidyr::crossing(freq_df) %>%
      mutate(Fraction = count / sum_of_cells)
  }
  
  plot <- data %>%
    ggplot(cell_type_fractions, aes(x = X, y = Fraction, fill = X)) +
    geom_boxplot() +
    facet_wrap(~ Split, scales = "free_y") +
    labs(
      title = "Fraction of Cell Types Within Each Sample and Tissue",
      x = "Group",
      y = "Fraction"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  
  print(plot)
  return (data)
}

ClusterComp3 <- function(
    seurat_object, 
    x = "Group", 
    hue = "orig.ident", 
    split = "cell_types",
    scales = "free", 
    reduce = F,
    return_splits = NA,
    remove_cells = 0
) {
  
  meta_df <- seurat_object@meta.data
  
  meta_df$X <- meta_df[[x]]
  meta_df$Hue  <- meta_df[[hue]]
  meta_df$Split <- meta_df[[split]]
  
  meta_df$Id <- meta_df[[hue]]
  if (reduce) {
  }

    id_counts <- meta_df %>% 
    dplyr::count(Id, name = "n_id")
  
  counts_df <- meta_df %>%
    dplyr::count(X, Id, Split, Hue, name = "count")
  
  data <- counts_df %>%
    left_join(id_counts, by = "Id") %>%
    group_by(Id) %>%
    mutate(new_count = ifelse(count > remove_cells, count, 0)) %>%
    mutate(count_diff = sum(count)-sum(new_count)) %>%
    mutate(n_id_2 = n_id - count_diff) %>%
    mutate(Fraction = new_count / n_id_2) %>%
    ungroup()
  
  plot <- data %>%
    filter(if (all(is.na(return_splits))) TRUE else Split %in% return_splits) %>%
    ggplot(aes(x = X, y = Fraction, fill = X)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1) + 
    facet_wrap(~ Split, scales = scales) +
    labs(
      title = "Fraction of Cell Types Within Each Sample (orig.ident)",
      x = x,
      y = "Fraction"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  
  print(plot)
  data <- list(
    "plot" = plot,
    "data" = data
  )
  return(data)
}

ClusterComp4 <- function(
    seurat_object, 
    x = "responder", 
    hue = "orig.ident", 
    split = "cell_types",
    pair = "status",
    scales = "free", 
    Fraction  = T,
    sample_id = "orig.ident",
    reduce = TRUE
) {
  
  meta_df <- seurat_object@meta.data
  
  # Keep the original columns for plotting
  meta_df$X <- meta_df[[x]]
  meta_df$Hue  <- meta_df[[hue]]
  meta_df$Split <- meta_df[[split]]
  meta_df$Pair <- meta_df[[pair]]
  
  # Use this column (Id) to compute fractions per sample
  meta_df$Id <- meta_df$Hue
  if (reduce) {
    # 1. Count total cells per Id (sample)
    total_counts <- meta_df %>% 
      dplyr::count(Id, name = "n")
    # 2. Determine the minimum cell count across Id
    min_n <- min(total_counts$n)
    
    # 3. Compute scaling factor per Id
    scaling_factors <- total_counts %>%
      mutate(scaling_factor = min_n / n) %>%
      dplyr::select(Id, scaling_factor)
    
    # 4. Count cells by (X, Id, Split, Hue)
    counts_df <- meta_df %>%
      dplyr::count(X, Id, Split, Hue, Pair, name = "count")
    
    # 5. Join scaling factors & compute fraction
    data <- counts_df %>%
      left_join(scaling_factors, by = "Id") %>%
      mutate(scaled_count = count * scaling_factor) %>%
      group_by(Id) %>%
      mutate(Fraction = scaled_count / sum(scaled_count)) %>%
      ungroup()
    
  } else {
    # Simple fraction without scaling
    # 1. Count total cells per Id
    id_counts <- meta_df %>% 
      count(Id, name = "n_id")
    
    # 2. Count cells by (X, Id, Split, Hue)
    counts_df <- meta_df %>%
      count(X, Id, Split, Hue, name = "count")
    
    # 3. Compute fraction of cells within each Id
    data <- counts_df %>%
      left_join(id_counts, by = "Id") %>%
      mutate(Fraction = count / n_id)
  }
  
  # Plot: x-axis is still X, but fraction is computed per Id
  if (Fraction){
    Y = "Fraction"
    label = "Proportion of Cells"
  }else{
    Y = "count"
    label = "Cell counts"
  }
  plot <- data %>%
    ggplot(aes(x = X,
               y = .data[[Y]],
               fill = Pair,
               group = interaction(X, Pair))) +
    # boxplots, dodged so baseline & week4 sit side by side
    geom_boxplot(position = position_dodge(width = 0.75), 
                 outlier.shape = NA) +
    # jittered points colored by patient (Hue), also dodged
    geom_jitter(position = position_jitterdodge(jitter.width = 0.15,
                                                dodge.width   = 0.75),
                size = 1,
                alpha = 0.7) +
    # one facet per cell type
    facet_wrap(~ Split, scales = "free_y") +
    scale_fill_brewer(palette = "Set2", name = "Timepoint") +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(
      title = label,
      x     = "Responder Status",
      y     = label
    ) +
    theme_minimal() +
    theme(
      axis.text.x  = element_text(angle = 30, hjust = 1),
      legend.position = "bottom"
    )
  
  print(plot)
  data <- list(
    "plot" = plot,
    "data" = data
  )
  return(data)
}
GO_KEGG <- function(seurat_object){
  
  
  deg.ls <- split(rownames(seurat_object), f = seurat_object$seurat_clusters)
  
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
    sample_name <- strsplit(file_paths[[i]],"/")[[1]][[3]]
    seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                     min.features = 100, 
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

ApplyDoubletFinder <- function (seurat_object){
  
  seurat_object$percent.mito <- PercentageFeatureSet(object = seurat_object, pattern = "^MT-")
  
  # run sctransform
  seurat_object <- SCTransform(seurat_object, vars.to.regress="percent.mito", assay="RNA",do.correct.umi=T,conserve.memory=T, verbose = T)
  
  seurat_object <- RunPCA(seurat_object)
  seurat_object <- RunUMAP(seurat_object, dims = 1:15)
  
  seurat_object <- FindNeighbors(seurat_object, dims = 1:15)
  
  seurat_object <- FindClusters(seurat_object, resolution = 0.3)
  
  sweep.res.list_seurat <- paramSweep(seurat_object, PCs = 1:15, sct = TRUE)
  sweep.stats_seurat <- summarizeSweep(sweep.res.list_seurat, GT = FALSE)
  bcmvn_seurat <- find.pK(sweep.stats_seurat)
  
  df <- data.frame(
    `cells.recovered` = c(800, 1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000),
    `multiplet.rate` = c(0.004, 0.008, 0.016, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076)
  )
  
  # 1. Fetch amount of cells from a seurat object
  num_cells <- length(seurat_object@meta.data$orig.ident)
  
  # 2. Check which row in the "# of Cells Recovered" column which this value is closest to
  closest_row_index <- which.min(abs(df$cells.recovered - num_cells))
  
  # 3. Define a new variable that matches the corresponding value on "Multiplet Rate (%)", that is on the same row.
  multiplet_rate <- df$multiplet.rate[closest_row_index]
  
  homotypic.prop <- modelHomotypic(seurat_object@meta.data$RNA_snn_res.0.3)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(multiplet_rate*nrow(seurat_object@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  seurat_object <- doubletFinder(seurat_object, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
  
  return (seurat_object)
}

CellCounts <- function(seurat_object, meta_data_var = "orig.ident", split = NULL, scales = "fixed") {
  # Extract the metadata from the Seurat object
  metadata <- seurat_object@meta.data
  
  if (is.null(split)) {
    # Count the number of cells per metadata var
    cell_counts <- metadata %>%
      group_by(.data[[meta_data_var]]) %>%
      summarise(n = n(), .groups = 'drop')
    
    print(cell_counts)
    
    # Calculate the overall mean cell count
    overall_mean <- mean(cell_counts$n, na.rm = TRUE)
    
    # Create the bar plot
    p <- ggplot(cell_counts, aes(x = .data[[meta_data_var]], y = n)) +
      geom_bar(stat = "identity") +
      theme_minimal() +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
      ) +
      labs(
        x = meta_data_var,
        y = "Cell Count",
        title = glue("Cell Count per {meta_data_var}")
      ) +
      geom_text(aes(label = n), vjust = -0.5, size = 3) +
      # Add overall mean line
      geom_hline(yintercept = overall_mean, color = "red", linetype = "dashed", size = 1) +
      # Add annotation for the mean
      annotate("text", 
               x = Inf, y = overall_mean, 
               label = paste("Mean =", round(overall_mean, 2)), 
               color = "red", 
               hjust = 1.1, vjust = -0.5, 
               size = 3.5)
    
  } else {
    # Check if split variable exists in metadata
    if (!split %in% colnames(metadata)) {
      stop(glue("Split variable '{split}' not found in the metadata."))
    }
    
    # Count the number of cells per metadata var and split
    cell_counts <- metadata %>%
      group_by(.data[[meta_data_var]], .data[[split]]) %>%
      summarise(n = n(), .groups = 'drop')
    
    print(cell_counts)
    
    # Calculate mean cell counts per split
    mean_per_split <- cell_counts %>%
      group_by(.data[[split]]) %>%
      summarise(mean_n = mean(n, na.rm = TRUE), .groups = 'drop')
    
    # Create the bar plot with facets
    p <- ggplot(cell_counts, aes(x = .data[[meta_data_var]], y = n)) +
      geom_bar(stat = "identity") +
      facet_wrap(as.formula(paste("~", split)), nrow = 3, scales = scales) +
      theme_minimal() +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
      ) +
      labs(
        x = meta_data_var,
        y = "Cell Count",
        title = glue("Cell Count per {meta_data_var} Split by {split}")
      ) +
      geom_text(aes(label = n), vjust = -0.5, size = 3) +
      # Add mean lines per split facet
      geom_hline(data = mean_per_split, aes(yintercept = mean_n), color = "red", linetype = "dashed", size = 1) +
      # Add annotations for the means
      geom_text(
        data = mean_per_split,
        aes(x = Inf, y = mean_n, label = paste("Mean =", round(mean_n, 2))),
        color = "red",
        hjust = 1.1,
        vjust = -0.5,
        size = 3.5
      )
    
  }
  
  return(p)
}

SCTNormalize <- function(seurat_object, vars.to.regress = NULL, n.features = 2000, vst.v = "v2"){
  
  # Get genes that arent Mitochondrial or chrY
  gene.table <- read.table("data/raw_external/hg38_gencode_v27.txt")
  non.sex.genes <- gene.table %>% 
    filter(!V2 %in% c("chrM", "chrY")) %>%
    pull(V1)
  print(vst.v)
  
  if (!is.null(vars.to.regress)){
    print(vars.to.regress)
    seurat_object <- SCTransform(seurat_object, do.correct.umi=T,conserve.memory=T, vars.to.regress = vars.to.regress, variable.features.n = n.features, vst.flavor = vst.v)
  }
  else(
    seurat_object <- SCTransform(seurat_object, do.correct.umi=T,conserve.memory=T, variable.features.n = n.features, vst.flavor = vst.v)
  )
  
  print("SCTransform done")
  gc()
  # Set 1s to 0s 
  print("Joining layers")
  #xr <- JoinLayers(seurat_object, assay = "RNA")[["RNA"]]$counts; gc()
  xr <- seurat_object[["RNA"]]$counts; gc()
  wr<-xr==0
  rm(xr);gc()
  print("Creating wr matrix")
  wr<-wr[rownames(seurat_object[["SCT"]]@counts),];gc()
  
  print("Creating xs matrix")
  xs<-as.matrix(seurat_object[["SCT"]]@counts);gc()
  
  print("Reset true 0s in xs")
  xs[as.matrix(wr)] <- 0
  #for (i in seq(1, ncol(xs), by = 1000)){
  #  xs[, i:min(i+9, ncol(xs))][as.matrix(wr[, i:min(i+9, ncol(xs))])] <- 0
  #  gc()
  #}
  rm(wr);gc()
  print("Reset true 0s in seurat object SCT counts")
  seurat_object[["SCT"]]@counts<-as.sparse(xs)
  rm(xs);gc()
  print("Overwrite SCT data slot in seurat object")
  #overwrite seurat_object slot by log2(seurat_object+1) instead of loge(seurat_object+1) from sctransform
  seurat_object[["SCT"]]@data<-as.sparse(log2(as.matrix(seurat_object[["SCT"]]@counts)+1))
  #genes <- c()
  #for (i in seq_along(seurat_object[['SCT']]@var.features)){
  #  if (seurat_object[['SCT']]@var.features[i] %in% non.sex.genes){
  #    len <- length(genes)
  #    genes[len + 1] <- seurat_object[['SCT']]@var.features[i]
  #  }
  #}
  #genes <- genes[genes %in% rownames(seurat_object)]
  #seurat_object[['SCT']]@var.features <- genes
  genes <- seurat_object[['SCT']]@var.features
  seurat_object <- RunPCA(seurat_object, features = genes)
  seurat_object <- RunUMAP(seurat_object, dims = 1:30)
  return(seurat_object)
}

CalcMetaGene <- function(seurat_object, gene_list){
  return(colMeans(x = seurat_object@assays$SCT$data[gene_list, ], na.rm = TRUE))
}

VizSummarizedLoadings <- function(seurat_object, dims = 15, genes = 20){
  loadings_matrix <- seurat_object@reductions$pca@feature.loadings[,1:dims]
  
  # Calculate absolute loadings across all PCs
  abs_loadings <- abs(loadings_matrix)
  
  # Sum up the absolute loadings for each feature
  aggregated_loadings <- rowSums(abs_loadings)
  
  # Rank features based on aggregated loadings
  sorted_features <- names(sort(aggregated_loadings, decreasing = TRUE))[1:genes]
  # Visualize aggregated loadings (example using ggplot2)
  
  df <- data.frame(Feature = sorted_features, Aggregated_Loadings = aggregated_loadings[sorted_features])
  df <- df[order(-df$Aggregated_Loadings), ]
  ggplot(df, aes(x = Feature, y = Aggregated_Loadings)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Features", y = "Aggregated Loadings", title = "Feature Importance Across PCs")
  
}

FilterCells <- function(seurat_object, gene_list, threshold = 1){
  for (i in seq_along(gene_list)){
    seurat_object <- seurat_object[, seurat_object[['RNA']]$counts[gene_list[i],] < threshold]
  }
  return(seurat_object)
}

CalculateAUCell <- function(seurat_object, gene_list, col.name){
  if (options()$Seurat.object.assay.version != 'v3'){
    seurat.option = options()$Seurat.object.assay.version
    options(Seurat.object.assay.version = 'v3')
  }else{
    seurat.option = 'v3'
  }
  gex_matrix <- seurat_object@assays$RNA$counts
  
  cells_rankings <- AUCell_buildRankings(gex_matrix, plotStats=FALSE, col.name)
  
  geneSets <- list(geneSet1=gene_list)
  
  cells_AUC <- AUCell_run(gex_matrix, geneSets, aucMaxRank=nrow(cells_rankings)*0.05)
  
  seurat_object <- AddMetaData(seurat_object, t(as.data.frame(cells_AUC@assays@data$AUC)), col.name = col.name)
  
  options(Seurat.object.assay.version = seurat.option)
  
  return(seurat_object)
}

AnnotateCell <- function(seurat_object, cell_names, group_ident = "cell_types", dimplot_ident = "orig.ident"){
  plot <- DimPlot(seurat_object, group.by = dimplot_ident)
  selected.cells <- CellSelector(plot, seurat_object)
  cells <- WhichCells(object = selected.cells, idents = "SelectedCells")
  cell.df <- data.frame(cells, cell_names, row.names = cells)
  cell.df$cells <- NULL
  seurat_object <- AddMetaData(seurat_object, cell.df, col.name = group_ident)
  return (seurat_object)
}

SubsetCells <- function(seurat_object, value = "orig.ident", n_cells = NULL){
  if (is.null(n_cells)){
    min_cells <- min(table(seurat_object[[value]]))
  }else if (n_cells == "median"){
    min_cells <- median(table(seurat_object[[value]]))
  }else if (is.numeric(n_cells)){
    min_cells <- n_cells
  }
  
  print(table(seurat_object[[value]]))
  
  seurat_object <- SetIdent(seurat_object, value = value)
  downsampled_cells <- lapply(unique(seurat_object[[value]])[[1]], function(ident) {
    cells <- colnames(seurat_object[,seurat_object[[value]] == ident])
    if (min_cells < length(cells)){
      n <- min_cells
    }else{
      n <- length(cells)
    }
    sample(cells, n)
  })
  # Flatten the list of downsampled cells into a single vector
  downsampled_cells <- unlist(downsampled_cells)
  # Subset the Seurat object
  seurat_object <- subset(seurat_object, cells = downsampled_cells)
  print(table(seurat_object[[value]]))
  return(seurat_object)
}



perform_enrichment <- function(genes, db = "GO_Biological_Process_2023", title = "", n_terms = 10) {
  Sys.sleep(1)
  df <- enrichr(genes, db)[[1]] %>%
    arrange(P.value) %>%
    slice_head(n = n_terms) %>%
    mutate(logP = -log10(P.value))
  plt <- df %>%
    ggplot(aes(x = reorder(Term, logP), y = logP)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(x = "", y = "-log10(P value)", title = paste("Top enriched GO terms ", title)) +
    theme_minimal()
  print(plt)
  return(df)
}

GO_Heatmap <- function(markers, n_terms = 5){
  
  # Convert gene names to ENTREZ IDs
  markers$entrez <- mapIds(org.Hs.eg.db, keys = markers$gene, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  
  
  # Split genes by cluster
  gene_list <- split(markers$entrez, markers$cluster)
  
  # Run GO enrichment for each cluster
  go_results <- lapply(gene_list, function(genes) {
    enrichGO(gene = genes[1:50], OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
             ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, maxGSSize = 50 )
  })
  # Convert GO results to a dataframe
  go_df <- do.call(rbind, lapply(names(go_results), function(cluster) {
    df <- as.data.frame(go_results[[cluster]])
    df$Cluster <- cluster
    return(df)
  }))
  
  # Select top GO terms per cluster
  top_go <- go_df %>%
    group_by(Cluster) %>%
    top_n(n_terms, wt = -log10(p.adjust)) %>% 
    ungroup()
  
  # Prepare matrix for heatmap (rows = GO terms, columns = clusters)
  go_matrix <- top_go %>%
    dplyr::select(Cluster, Description, p.adjust) %>%
    mutate(log_pval = -log10(p.adjust)) %>%
    tidyr::spread(Cluster, log_pval, fill = 0) %>%
    dplyr::select(-p.adjust)
  
  # Convert to matrix
  heatmap_matrix <- as.matrix(go_matrix[,-1]) 
  rownames(heatmap_matrix) <- go_matrix$Description
  
  ph <- pheatmap(heatmap_matrix,
                 cluster_rows = TRUE, 
                 cluster_cols = TRUE, 
                 color = colorRampPalette(c("white","red"))(5),
                 fontsize_row = 8, name = "-log10(p)")
  print(ph)
  return(top_go)
}

CorrelateMarkers <- function(seurat_obj, markerdf, n_filt = NULL, seurat = T){
  if (!is.null(n_filt)){
    markerdf<-markerdf[c(1:30),]
  }
  
  
  #seurat object in which we estimate all correlations
  #data<-readRDS("/Users/jonas_sjolund/Documents/Jonas/Lund Research/Datasets/CAFs_Cords_Bodenmiller/scRNA-seq/PDAC_gse21_fibro_tumour.rds")
  data <- seurat_obj
  
  
  ##########################################################################################
  #Define functions
  ##########################################################################################
  getZscore <- function(dat){
    z.dat <- matrix(nrow=nrow(dat), ncol=ncol(dat), dimnames=list(rownames(dat),colnames(dat)))
    
    dat.mean <- apply(dat, 1, mean)
    dat.stdev <- apply(dat, 1, sd)
    
    for(i in 1:nrow(dat)){
      x <- as.numeric(dat[i,])
      z.dat[i,] <- (x - dat.mean[i])/dat.stdev[i]
    }
    z.dat <- as.data.frame(z.dat)
    return(z.dat)
  }
  
  computeScore = function(sign){
    if (seurat == T){
      sign.dat <- GetAssayData(object = data,slot ="data")
    }else{
      sign.dat <- data
    }
    print(rownames(sign.dat))
    sign.dat <- sign.dat[which(rownames(sign.dat)%in%sign),]
    
    sign.dat <- as.matrix(sign.dat)
    print(dim(sign.dat))
    print(names(sign))
    
    # Z-score
    z.score.dat <- getZscore(sign.dat)
    z.score.avg <- apply(z.score.dat,2,mean,na.rm=TRUE)
  }
  
  
  ##########################################################################################
  ##########################################################################################
  
  # Apply get_score to each column and combine results as columns
  result_list <- map(markerdf, computeScore)  # Apply function to each column
  score_all <- bind_cols(result_list)  # Combine results column-wise
  #correlation matrix
  method="pearson"
  mat.cor = cor(score_all,method=method, use="pairwise.complete.obs")
  
  col2 = colorRampPalette(c("#046380","white","#B9121B"))
  col = col2(200)
  library(corrplot)
  dev.off()
  corrplot(mat.cor,
           order="hclust",
           hclust.method="ward.D2",
           tl.cex=0.65,
           tl.col="black",
           col=col,method = "square",)
  
}

MultipleModuleScores <- function(seurat_obj, markerdf, mode = "Seurat"){
  if (mode != "Seurat"){
    for (column in colnames(markerdf)){
      print(column)
      seurat_obj <- AddModuleScore(seurat_obj, list(markerdf[[column]]), name = column)
    }}
  else{
    for (name in unique(markerdf$cluster)){
      genes <- markerdf[markerdf$cluster == name,]$gene
      seurat_obj <- AddModuleScore(seurat_obj, list(genes), name = name)
    }
  }
  return(seurat_obj)
}





GO_Heatmap_EnrichR <- function(markers, n_terms = 5, n_genes = 50, up = T, max.cutoff = 5){
  # Define the database to use for enrichment (GO Biological Process)
  db <- "GO_Biological_Process_2023"
  
  if (up){
    markers <- markers %>% dplyr::group_by(gene) %>% 
      arrange(desc(avg_log2FC)) %>% 
      slice_head(n = 1) %>% 
      ungroup()
    heat_col <- "red"
  } else {
    markers <- markers %>% dplyr::group_by(gene) %>% 
      arrange(avg_log2FC) %>% 
      slice_head(n = 1) %>% 
      ungroup()
    heat_col <- "blue"
  }
  
  # Split genes by cluster
  gene_list <- split(markers$gene, markers$cluster)
  
  # Run GO enrichment for each cluster using enrichR
  go_results <- lapply(names(gene_list), function(cluster) {
    print(cluster)
    print(gene_list[[cluster]][1:n_genes])
    enriched <- enrichr(gene_list[[cluster]][1:n_genes], db)[[db]]
    enriched$Cluster <- cluster  # Add cluster identifier
    Sys.sleep(1)
    return(enriched)
  })
  
  # Convert GO results to a dataframe
  go_df <- bind_rows(go_results)
  
  # Select top GO terms per cluster
  top_go <- go_df %>%
    group_by(Cluster) %>%
    slice_max(order_by = -log10(Adjusted.P.value), n = n_terms, with_ties = FALSE) %>%
    ungroup()
  
  # Prepare matrix for heatmap (rows = GO terms, columns = clusters)
  go_matrix <- top_go %>%
    dplyr::select(Cluster, Term, Adjusted.P.value) %>%
    mutate(log_pval = -log10(Adjusted.P.value)) %>%
    tidyr::spread(Cluster, log_pval, fill = 0) %>%
    dplyr::select(-Adjusted.P.value)
  
  heatmap_matrix <- as.matrix(go_matrix[,-1])
  rownames(heatmap_matrix) <- go_matrix$Term
  # Convert matrix to dataframe for manipulation
  heatmap_df <- as.data.frame(heatmap_matrix)
  
  # Add row names as a column for grouping
  heatmap_df$GO_Term <- rownames(heatmap_matrix)
  
  # Sum values for duplicate GO terms
  heatmap_df <- heatmap_df %>%
    group_by(GO_Term) %>%
    summarise(across(everything(), sum, na.rm = TRUE))
  
  # Convert back to matrix
  heatmap_matrix <- as.matrix(heatmap_df[,-1])
  rownames(heatmap_matrix) <- heatmap_df$GO_Term
  heatmap_numbers <- round(heatmap_matrix, 2)
  
  # Define the color palette
  colors <- colorRampPalette(c("white", heat_col))(50)
  
  # Set breaks to range from 0 to 5
  breaks <- seq(0, max.cutoff, length.out = length(colors) + 1)
  
  # Plot
  ph <- pheatmap(heatmap_matrix,
                 cluster_rows = TRUE, 
                 cluster_cols = TRUE, 
                 color = colors,
                 breaks = breaks,
                 fontsize_row = 8,
                 name = "-log10(p)",
                 display_numbers = heatmap_numbers)
  
  print(ph)
  return(top_go)
}

save_last_plot <- function(fig_name,
                           out_dir,
                           width  = 7,   
                           height = 7) {  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  filepath <- file.path(out_dir, paste0(fig_name, ".pdf"))
  dev.copy(pdf,
           file   = filepath,
           width  = width,
           height = height)
  dev.off()
  message("Saved last plot to: ", filepath,
          sprintf("  [%0.1f×%0.1f inches]", width, height))
}


plotSpecificIdent <- function(
    object,
    meta.column,
    ident.level,
    reduction   = "umap",
    dims        = c(1,2),
    cols        = c("lightgrey","red"),
    pt.size     = 0.5,
    ...
) {
  # 1) Pull out metadata and set as active identities
  if ( ! meta.column %in% colnames(object@meta.data) ) {
    stop("‘", meta.column, "’ not found in object@meta.data")
  }
  # Create a factor of length cells that is either “other” or the target level
  meta.vec <- object@meta.data[[meta.column]]
  highlight <- ifelse(meta.vec == ident.level, as.character(ident.level), "other")
  highlight <- factor(highlight, levels = c("other", as.character(ident.level)))
  
  # 2) Run DimPlot but override grouping to use our two-level factor
  p <- Seurat::DimPlot(
    object,
    reduction = reduction,
    group.by  = NULL,           # ignore default identities
    cells     = Cells(object),  # include all cells
    pt.size   = pt.size,
    # supply our own labels & colors
    cols      = cols,
    ...
  )
  
  # 3) Since DimPlot uses the active identities internally, we can stash & swap
  original.idents <- Seurat::Idents(object)
  Seurat::Idents(object) <- highlight
  p <- Seurat::DimPlot(
    object,
    reduction = reduction,
    group.by  = NULL,
    cells     = Cells(object),
    pt.size   = pt.size,
    cols      = cols,
    dims      = dims,
    ...
  ) + ggplot2::labs(title = paste0(meta.column, " == ", ident.level))
  # restore
  Seurat::Idents(object) <- original.idents
  return(p)
}

DotPlot3 <- function(){
  df <- data.frame(
    category   = factor(c("A", "B", "C", "D", "E"),
                        levels = c("A","B","C","D","E")),
    value      = c(10, 23, 17, 9, 28),
    size_var   = c(50, 200, 100, 80, 300),
    colour_var = c(5.1, 3.3, 7.8, 2.2, 9.5)   # some continuous measure
  )
  
  ggplot(df, aes(
    x     = value,
    y     = category,
    size  = size_var,
    colour= colour_var     # map your continuous variable here
  )) +
    geom_point(alpha = 0.8) +
    # make area proportional to `size_var`
    scale_size_area(max_size = 15) +
    # choose a continuous colour gradient:
    scale_colour_gradient(
      low  = "lightblue",
      high = "darkblue",
      name = "Colour\nVariable"
    ) +
    # or use a perceptually uniform Viridis palette:
    # scale_colour_viridis_c(name = "Colour\nVariable") +
    labs(
      title = "Bubble Plot: Size & Colour by Two Continuous Variables",
      x     = "Value",
      y     = "Category"
    ) +
    theme_minimal(base_size = 14)
}

FilterGenes <- function(seurat_object){
  genes.use <- grep(pattern = "^RP[SL][[:digit:]]|^RP[[:digit:]]|^RPSA|^MT-|MRP[SL]",
                    rownames(seurat_object),
                    value=TRUE, invert=TRUE)
  DefaultAssay(seurat_object ) <- "RNA"
  seurat_object@assays$SCT <- NULL
  seurat_object <- seurat_object[genes.use,]
  
}



CorrelateCells <- function(df, 
                           nclust = 1, 
                           seurat_object = NA, 
                           features = NA, 
                           split = "seurat_clusters",
                           proportions = NA
                           ){
  # --- base matrix for the correlation heatmap ---
  matrix <- df %>% 
    dplyr::select(Split, Fraction, Id) %>%
    tidyr::pivot_wider(names_from = Split, values_from = Fraction, values_fill = 0) %>%
    as.matrix()
  sample_vector <- matrix[,1]
  matrix <- matrix[,-1]
  matrix <- as.matrix(apply(matrix, 2, as.numeric))
  rownames(matrix) <- sample_vector
  
  # correlation matrix (rows/cols = Splits)
  corr_mat <- cor(matrix, method = "spearman")
  
  id_levels <- if (is.factor(df$Id)) levels(df$Id) else sort(unique(df$Id))
  
  id_comp_df <- df %>%
    dplyr::group_by(Split, Id) %>%
    dplyr::summarise(frac = sum(Fraction, na.rm = TRUE), .groups = "drop") %>%
    group_by(Split) %>%
    tidyr::nest() %>%
    dplyr::mutate(vals = purrr::map(data, ~setNames(.x$frac, .x$Id))) %>%
    dplyr::pull(vals)
  
  
  frame_bar_anno <- function(index) {
    rn <- rownames(corr_mat)[index]
    v <- prop_aligned[rn]
    v[!is.finite(v)] <- 0
    v <- pmin(pmax(v, 0), 1)
    n <- length(v)
    if (n == 0) return(invisible())
    h <- 1/n
    for (i in seq_len(n)) {
      y0 <- (n - i) / n
      grid::grid.rect(x = 0, y = y0, width = 1, height = h,
                      just = c("left", "bottom"),
                      gp = grid::gpar(fill = NA, col = "black", lwd = 0.6))
      if (v[i] > 0) {
        grid::grid.rect(x = 0, y = y0, width = v[i], height = h,
                        just = c("left", "bottom"),
                        gp = grid::gpar(fill = "grey50", col = NA))
      }
    }
  }

    id_levels <- if (is.factor(df$Id)) levels(df$Id) else sort(unique(df$Id))
  
  id_comp_mat <- df %>%
    dplyr::group_by(Split, Id) %>%
    dplyr::summarise(frac = sum(Fraction, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Id, values_from = frac, values_fill = 0) %>%
    tibble::column_to_rownames("Split") %>%
    as.matrix()
  
  row_sums <- rowSums(id_comp_mat)
  row_sums[row_sums == 0] <- 1
  id_comp_mat <- id_comp_mat / row_sums
  
  id_comp_mat <- id_comp_mat[rownames(corr_mat), id_levels, drop = FALSE]
  
  id_cols <- setNames(grDevices::hcl.colors(length(id_levels), "Zissou1"), id_levels)
  
  ra_left <- ComplexHeatmap::rowAnnotation(
    "Sample composition" = ComplexHeatmap::anno_barplot(
      id_comp_mat,
      which = "row",
      border = FALSE,
      bar_width = 1,
      gp = grid::gpar(fill = id_cols[colnames(id_comp_mat)])
    ),
    annotation_name_rot = 90
  )
  
  ht <- ComplexHeatmap::Heatmap(
    corr_mat, name = "Spearman",
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = "euclidean",
    clustering_method_columns = "ward.D",
    clustering_method_rows = "ward.D",
    border = TRUE,
    #row_split = nclust,
    row_gap = grid::unit(5, "mm"),
    left_annotation = ra_left,
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid::grid.rect(x = x, y = y, width = width, height = height,
                      gp = grid::gpar(col = "black", fill = NA, lwd = 0.5))
    }
  )
  
  ComplexHeatmap::draw(ht)
}
top_correlated_genes <- function(seu, gene, assay = NULL, slot = "data",
                                 method = "pearson", n = 50, exclude_target = TRUE) {
  assay <- assay %||% DefaultAssay(seu)
  mat   <- GetAssayData(seu, assay = assay, slot = slot)
  
  if (!gene %in% rownames(mat)) stop("Gene not found in assay")
  
  # pull target vector (numeric, dense)
  y <- as.numeric(mat[gene, ])
  
  # transpose mat so genes are columns, cells are rows
  cors <- cor(y, as.matrix(t(mat)), method = method, use = "pairwise.complete.obs")
  
  res <- data.frame(
    gene = rownames(mat),
    cor  = as.numeric(cors),
    stringsAsFactors = FALSE
  )
  
  if (exclude_target) res <- res[res$gene != gene, ]
  
  res[order(res$cor, decreasing = TRUE), ][1:n, ]
}




seurat_to_spe <- function(seu, sample_id, img_id) {
  ## Convert to SCE
  sce <- Seurat::as.SingleCellExperiment(seu)
  
  ## Extract spatial coordinates
  spatialCoords <- as.matrix(
    seu@images[[img_id]]@coordinates[, c("imagecol", "imagerow")])
  
  ## Extract and process image data
  img <- SpatialExperiment::SpatialImage(
    x = as.raster(seu@images[[img_id]]@image))
  
  imgData <- DataFrame(
    sample_id = sample_id,
    image_id = img_id,
    data = I(list(img)),
    scaleFactor = seu@images[[img_id]]@scale.factors$lowres)
  
  # Convert to SpatialExperiment
  spe <- SpatialExperiment(
    assays = assays(sce),
    rowData = rowData(sce),
    colData = colData(sce),
    metadata = metadata(sce),
    reducedDims = reducedDims(sce),
    altExps = altExps(sce),
    sample_id = sample_id,
    spatialCoords = spatialCoords,
    imgData = imgData
  )
  # indicate all spots are on the tissue
  spe$in_tissue <- 1
  spe$sample_id <- sample_id
  # Return Spatial Experiment object
  spe
}


HighlightCell <- function(
    seu_obj,
    cell_ident = NULL,
    threshold = 0,
    feature = NULL,
    group = NULL,
    fov = NULL,
    plot_type = "Dim",
    assay = NULL,          
    slot = "data",
    new_group = NULL
) {
  if (is.null(assay)) assay <- DefaultAssay(seu_obj)
  
  if (!is.null(feature)) {
    if (sum(feature %in% rownames(seu_obj)) == length(feature)) {
      vals <- as.numeric(GetAssayData(seu_obj, assay = assay, slot = slot)[feature, , drop = TRUE])
      names(vals) <- colnames(seu_obj)
    } else if (feature %in% colnames(seu_obj@meta.data)) {
      vals <- seu_obj@meta.data[[feature]]
      names(vals) <- rownames(seu_obj@meta.data)
    } else {
      stop(sprintf("`%s` not found as a gene in assay '%s' (slot '%s') or as a metadata column.",
                   feature, assay, slot))
    }
    cells <- names(vals)[vals > threshold]
  } else if (!is.null(group) && !is.null(cell_ident)) {
    if (!group %in% colnames(seu_obj@meta.data)) {
      stop(sprintf("Group column `%s` not found in meta.data.", group))
    }
    cells <- unique(rownames(seu_obj@meta.data)[seu_obj@meta.data[[group]] == cell_ident])
  } else {
    stop("Provide either `feature` (with optional `threshold`) or both `group` and `cell_ident`.")
  }
  if (length(cells) == 0) {
    message("No cells passed the filter; nothing to highlight.")
  } else {
    message(sprintf("Highlighting %d cells.", length(cells)))
  }
  
  
  if (plot_type == "Dim") {
      print(DimPlot(seu_obj, 
              cells.highlight = cells
              ))
  } else if (plot_type == "Spatial") {
    print(ImageDimPlot(seu_obj,
                 cells = cells, 
                 fov = fov,
                 group.by = "cell_types_2"))
  } else if (plot_type == "Feature") {
      FeaturePlot(seu_obj, 
                  features = feature, 
                  cells = colnames(seu_obj)
                  )
  } else {
      stop("plot_type must be one of 'Dim', 'Spatial', or 'Feature'.")
  }
  return(cells)
}

drawInteraction <- function(interaction_table, title =""){
  interaction_df <- as.data.frame(interaction_table)
  
  # Keep only non-zero interactions
  interaction_df <- interaction_df[interaction_df$Freq > 0, ]
  
  # 1) Give the canvas some breathing room + silence overflow warnings
  circos.clear()
  circos.par(
    canvas.xlim = c(-1.2, 1.2),
    canvas.ylim = c(-1.5, 1.5),
    points.overflow.warning = FALSE,
    track.margin = c(0.01, 0.01),
    gap.degree = 4
  )
  
  # 2) Draw the chord diagram but DON'T draw default names
  chordDiagram(
    interaction_df,
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.16) # taller label track
  )
  
  # 3) Add your own rotated labels, nudged inward so they don't clip
  circos.trackPlotRegion(
    track.index = 1,
    panel.fun = function(x, y) {
      sector = get.cell.meta.data("sector.index")
      xcenter = get.cell.meta.data("xcenter")
      ylim    = get.cell.meta.data("ylim")
      circos.text(
        x = xcenter,
        y = ylim[1] + mm_y(2.5),   # nudge inward; increase if still clipping
        labels = sector,
        facing = "clockwise",      # or "reverse.clockwise", "inside"
        niceFacing = TRUE,
        adj = c(0, 0.5),
        cex = 0.8
      )
    },
    bg.border = NA
  )
  title(title, cex.main = 1.2)
  
}

ProsegToSeurat <- function(proseg_output_path, expected_counts_basename="expected-counts", cell_metadata_basename="cell-metadata") {
  expected_counts_path <- file.path(proseg_output_path, paste0(expected_counts_basename, ".csv.gz"))
  if (!file.exists(expected_counts_path)) {
    expected_counts_path <- file.path(proseg_output_path, paste0(expected_counts_basename, ".parquet"))
    if (!file.exists(expected_counts_path)) {
      stop("Can't find expected-counts file.")
    }
    expected_counts <- arrow::read_parquet(expected_counts_path)
    
  } else {
    expected_counts <- readMM(gzfile(expected_counts_path))
  }
  
  cell_metadata_path <- file.path(proseg_output_path, paste0(cell_metadata_basename, ".csv.gz"))
  if (!file.exists(cell_metadata_path)) {
    cell_metadata_path <- file.path(proseg_output_path, paste0(cell_metadata_basename, ".parquet"))
    if (!file.exists(cell_metadata_path)) {
      stop("Can't find cell-metadata file.")
    }
    cell_metadata <- arrow::read_parquet(cell_metadata_path)
    
  } else {
    cell_metadata <- read.csv(cell_metadata_path, header=TRUE, sep=",")
  }
  
  # exclude cells with no assigned transcripts that end up with undefinied centroids
  mask <- is.finite(cell_metadata$centroid_x) & is.finite(cell_metadata$centroid_y)
  cell_metadata <- cell_metadata[mask,]
  expected_counts <- expected_counts[mask,]
  
  sobj <- Seurat::CreateSeuratObject(
    counts=Matrix::Matrix(t(as.matrix(expected_counts)), sparse=TRUE),
    meta.data=cell_metadata,
    assay="RNA"
  )
  
  # From: https://github.com/satijalab/seurat/issues/2790#issuecomment-721400652
  coords_df <- cell_metadata[c("centroid_x", "centroid_y")]
  names(coords_df) <- c("x", "y")
  rownames(coords_df) <- colnames(sobj)
  
  sobj@images$image =  new(
    Class = 'SlideSeq',
    assay = "Spatial",
    key = "image_",
    coordinates =  coords_df)
  
  return(sobj)
}

GetSampleFov <- function(seu, fov) {
  cells_in_fov <- Cells(seu[[fov]])
  sample_id <- unique(seu$sample[cells_in_fov])
  print(paste("Sample identity for cells in", fov, "is:", sample_id))
  return(sample_id)
}

CellDensity <- function(seu, ct1, fov, ct2 = NA, rotate = T, mirror = F){
  sample_name <- GetSampleFov(seu, fov)
  p <- ImageDimPlot(seu, fov = fov, alpha = 0.3, group.by = "cell_types", border.size = 0)
  coords <- GetTissueCoordinates(seu, which = "centroids", image = fov)
  coords$cell_type <- seu[,seu$sample == sample_name]$cell_types

  macro_coords <- coords[coords$cell_type == ct1, ]
  macro_coords2 <- coords[coords$cell_type == ct2, ]
  
  if (rotate){
    x <- macro_coords$x
    y <- macro_coords$y
    x2 <- macro_coords2$x
    y2 <- macro_coords2$y
    
    macro_coords$y <- x
    macro_coords$x <- y
    macro_coords2$y <- x2
    macro_coords2$x <- y2
  }  
  if (mirror){
    x <- macro_coords$x
    y <- macro_coords$y
    x2 <- macro_coords2$x
    y2 <- macro_coords2$y
    
    macro_coords$y <- x
    macro_coords$x <- y
    macro_coords2$y <- x2
    macro_coords2$x <- y2
  }  
  
  
  if (is.na(ct2)){
    p + 
      new_scale_fill() + 
      stat_density_2d(
        data = macro_coords, 
        aes(x = x, y = y, fill = after_stat(level)), 
        geom = "polygon", 
        alpha = 0.8, 
        bins = 10,
        h = c(100, 100),
        inherit.aes = FALSE
      ) +
      scale_fill_viridis_c(
        option = "magma",
        name = paste(ct1, "\nDensity")
      )
  }else{
    p + 
      new_scale_fill() + 
      stat_density_2d(
        data = macro_coords, 
        aes(x = x, y = y, fill = after_stat(level)), 
        geom = "polygon", 
        alpha = 0.8, 
        bins = 10,
        h = c(100, 100),
        inherit.aes = FALSE
      ) +
      scale_fill_viridis_c(
        option = "magma",
        name = paste(ct1, "\nDensity")
      ) +
      new_scale_fill() + 
      stat_density_2d(
        data = macro_coords2, 
        aes(x = x, y = y, fill = after_stat(level)), 
        geom = "polygon", 
        alpha = 0.8, 
        bins = 10,
        h = c(100, 100),
        inherit.aes = FALSE
      ) +
      scale_fill_viridis_c(
        option = "cividis",
        name = paste(ct2, "\nDensity")
      )
  }
}

CalculateNeighbourhoods <- function(xenium.obj, k = 10, nclust = 10){
  
  prop_list <- list()
  for (fov in Images(xenium.obj)){
    message(paste("Processing fov: ",fov))
    df_cells <- GetTissueCoordinates(xenium.obj, image = fov)
    celltypes <- data.frame(xenium.obj$cell_types)
    celltypes$cell <- rownames(celltypes)
    df_cells <- df_cells %>% left_join(y = celltypes, by = "cell")
    coords <- df_cells %>%
      dplyr::select(cell, x, y, xenium.obj.cell_types)
    
    coord_mat <- coords %>%
      dplyr::select(x, y) %>%
      as.matrix()
    
    rownames(coord_mat) <- coords$cell
    
    
    
    knn <- get.knn(coord_mat, k = k)
    celltypes <- coords$xenium.obj.cell_types
    names(celltypes) <- coords$cell
    
    neigh_prop_list <- lapply(seq_len(nrow(knn$nn.index)), function(i) {
      neigh_cells <- rownames(coord_mat)[knn$nn.index[i, ]]
      neigh_types <- celltypes[neigh_cells]
      prop.table(table(neigh_types))
    })
    
    all_types <- sort(unique(unlist(lapply(neigh_prop_list, names))))
    
    neigh_prop_mat <- matrix(
      0,
      nrow = length(neigh_prop_list),
      ncol = length(all_types),
      dimnames = list(names(neigh_prop_list) %||% NULL, all_types)
    )
    
    for (i in seq_along(neigh_prop_list)) {
      v <- neigh_prop_list[[i]]
      neigh_prop_mat[i, names(v)] <- as.numeric(v)
    }
    
    rownames(neigh_prop_mat) <- coords$cell
    prop_list[[fov]] <- neigh_prop_mat
    gc()
  }
  all_cols <- Reduce(union, lapply(prop_list, colnames))
  
  
  # pad each matrix with missing columns
  prop_list_filled <- lapply(prop_list, function(m) {
    missing <- setdiff(all_cols, colnames(m))
    if (length(missing) > 0) {
      m <- cbind(m, matrix(0, nrow(m), length(missing),
                           dimnames = list(rownames(m), missing)))
    }
    m[, all_cols, drop = FALSE]
  })
  
  prop_mat <- do.call(rbind, prop_list_filled)
  
  
  kmeans_model <- MiniBatchKmeans(prop_mat, 
                                  clusters = nclust,
                                  batch_size = 50000,
                                  num_init = 5,
                                  max_iters = 100,
                                  early_stop_iter = 10,
                                  initializer = "kmeans++",
                                  seed = 7991) 
  cluster_assignments <- predict_KMeans(prop_mat, kmeans_model$centroids)  
  #prop_mat$cluster <- as.factor(cluster_assignments)
  glimpse(prop_mat)
  centroids <- kmeans_model$centroids
  rownames(centroids) <- paste0("K", seq_len(nrow(centroids)))
  colnames(centroids) <- colnames(prop_mat)
  centroids <- t(t(centroids) / colMaxs(centroids))
  ht <- Heatmap(
    centroids,
    name = "prop",
    row_title = "K-means clusters",
    column_title = "Neighbor cell-type proportions",
    cluster_rows = F,
    col = colorRamp2(c(0,1),
                     c("white", "red")
    ),
    cluster_columns = TRUE
  )
  
  cluster_assignments <- predict_KMeans(prop_mat, kmeans_model$centroids)
  
  cluster_assignments_named <- setNames(
    cluster_assignments,
    rownames(prop_mat)
  )
  
  xenium.obj <- AddMetaData(
    object = xenium.obj,
    metadata = cluster_assignments_named,
    col.name = "kmeans_cluster"
  )
  return(list("xenium" = xenium.obj, "ht" = ht))
}

choose_k_proportional_stability <- function(coord_mat, celltypes, k_grid = 5:80, step = 5) {
  k_grid <- k_grid[k_grid %% step == 0]
  mats <- lapply(k_grid, \(k) neighbor_prop_mat(coord_mat, celltypes, k))
  
  changes <- sapply(seq_along(mats)[-1], function(i){
    mean(rowSums(abs(mats[[i]] - mats[[i-1]])))
  })
  
  df <- data.frame(
    k = k_grid[-1],
    mean_L1_change = changes
  )
  
  
  thr <- quantile(df$mean_L1_change, 0.2)  # heuristic: lower 20% of changes
  k_star <- df$k[min(which(df$mean_L1_change <= thr))]
  
  list(k_star = k_star, curve = df)
}

choose_k_predictive <- function(coord_mat, celltypes, k_grid = 5:80) {
  y <- factor(celltypes[rownames(coord_mat)])
  ks <- k_grid
  
  acc <- sapply(ks, function(k){
    P <- neighbor_prop_mat(coord_mat, y, k)   # uses neighbor labels only
    pred <- colnames(P)[max.col(P, ties.method = "first")]
    mean(pred == as.character(y))
  })
  
  df <- data.frame(k = ks, accuracy = acc)
  k_star <- df$k[which.max(df$accuracy)]
  
  list(k_star = k_star, curve = df)
}

choose_k_predictive_spatial_cv <- function(coord_mat,
                                           celltypes,
                                           k_grid = 5:80,
                                           n_folds = 5,
                                           seed = 1,
                                           n_blocks_per_axis = 6) {
  if (!requireNamespace("RANN", quietly = TRUE)) {
    stop("Please install RANN: install.packages('RANN')")
  }
  
  y <- factor(celltypes[rownames(coord_mat)])
  if (any(is.na(y))) stop("Some coord_mat rownames are missing in `celltypes`.")
  if (ncol(coord_mat) < 2) stop("coord_mat should have at least 2 columns (x,y).")
  
  set.seed(seed)
  n <- nrow(coord_mat)
  classes <- levels(y)
  
  xy <- coord_mat[, 1:2, drop = FALSE]
  # Create block ids by discretizing x and y into n_blocks_per_axis bins each
  bx <- cut(xy[, 1], breaks = n_blocks_per_axis, labels = FALSE, include.lowest = TRUE)
  by <- cut(xy[, 2], breaks = n_blocks_per_axis, labels = FALSE, include.lowest = TRUE)
  block_id <- interaction(bx, by, drop = TRUE)
  
  # Assign blocks (not cells) to folds
  blocks <- levels(block_id)
  block_fold <- sample(rep(seq_len(n_folds), length.out = length(blocks)))
  names(block_fold) <- blocks
  folds <- block_fold[as.character(block_id)]
  
  predict_with_k <- function(k, train_idx, test_idx) {
    train_xy <- coord_mat[train_idx, , drop = FALSE]
    test_xy  <- coord_mat[test_idx,  , drop = FALSE]
    train_y  <- y[train_idx]
    
    nn <- RANN::nn2(data = train_xy, query = test_xy, k = k)
    nn_idx <- nn$nn.idx
    neigh_lab <- matrix(train_y[nn_idx], nrow = nrow(nn_idx), ncol = ncol(nn_idx))
    
    pred <- apply(neigh_lab, 1, function(labs) {
      counts <- tabulate(match(labs, classes), nbins = length(classes))
      classes[which.max(counts)]
    })
    
    pred
  }
  
  acc <- sapply(k_grid, function(k) {
    fold_acc <- numeric(n_folds)
    for (f in seq_len(n_folds)) {
      test_idx  <- which(folds == f)
      train_idx <- which(folds != f)
      
      # If a fold is empty (can happen with few blocks), skip safely
      if (length(test_idx) == 0 || length(train_idx) == 0) {
        fold_acc[f] <- NA_real_
        next
      }
      
      pred <- predict_with_k(k, train_idx, test_idx)
      fold_acc[f] <- mean(pred == as.character(y[test_idx]))
    }
    mean(fold_acc, na.rm = TRUE)
  })
  
  df <- data.frame(k = k_grid, spatial_cv_accuracy = acc)
  k_star <- df$k[which.max(df$spatial_cv_accuracy)]
  
  list(k_star = k_star, curve = df, folds = folds, block_id = block_id)
}


neighbor_prop_mat <- function(coord_mat, celltypes, k) {
  knn <- get.knn(coord_mat, k = k)
  all_types <- sort(unique(celltypes))
  
  out <- matrix(0, nrow(coord_mat), length(all_types),
                dimnames = list(rownames(coord_mat), all_types))
  
  for (i in seq_len(nrow(coord_mat))) {
    neigh_cells <- rownames(coord_mat)[knn$nn.index[i,]]
    tab <- table(celltypes[neigh_cells])
    out[i, names(tab)] <- as.numeric(tab) / sum(tab)
  }
  out
}


concatenate_fovs = function(obj,
                            assay = "RNA",
                            n_cols = 2, 
                            offset = 5000,
                            fov_name = "combined", 
                            append = TRUE){
  if(!("Seurat" %in% .packages()))
  {
    library(Seurat)
  }
  
  #get shape of final image
  n_fovs = length(obj@images)
  n_rows = ceiling(n_fovs / n_cols)
  
  #make new coordinates for combined object
  starting_x = 0 #where to place leftmost point
  starting_y = 0 #where to place topmost point
  final_centroids = NULL #concatenated version
  final_molecules = NULL
  
  for(image in 1:length(obj@images))
  {
    #update centroid x and y based on previously defined offsets
    if(image == 1)
    {
      #uses metadata (nsides, theta, radius) from first fov only
      final_centroids = obj@images[[image]]$centroids
      updated_centroids = final_centroids
    } else
    {
      updated_centroids = obj@images[[image]]$centroids
      updated_centroids@coords[,"x"] = updated_centroids@coords[,"x"] + starting_x
      updated_centroids@coords[,"y"] = updated_centroids@coords[,"y"] + starting_y
      final_centroids@coords = rbind(final_centroids@coords, updated_centroids@coords)
      final_centroids@cells = c(final_centroids@cells, updated_centroids@cells)
      #expand out limits 
      #note that min is inherited from image 1
      final_centroids@bbox["x", "max"] = max(final_centroids@coords[, "x"])
      final_centroids@bbox["y", "max"] = max(final_centroids@coords[, "y"])
    }
    
    #update molecule locations
    if(image == 1)
    {
      final_molecules = obj@images[[image]]$molecules
    } else
    {
      updated_molecules = obj@images[[image]]$molecules
      for(molecule in names(updated_molecules)){
        updated_molecules[[molecule]]@coords[,"x"] = updated_molecules[[molecule]]@coords[,"x"] + starting_x
        updated_molecules[[molecule]]@coords[,"y"] = updated_molecules[[molecule]]@coords[,"y"] + starting_y
        updated_molecules[[molecule]]@bbox["x",] = updated_molecules[[molecule]]@bbox["x",] + starting_x
        updated_molecules[[molecule]]@bbox["y",] = updated_molecules[[molecule]]@bbox["y",] + starting_y
      }
      
      for(molecule in names(final_molecules)){
        final_molecules[[molecule]]@coords = rbind(final_molecules[[molecule]]@coords, updated_molecules[[molecule]]@coords)
        #update boundaries
        #note that min is inherited from image 1
        final_molecules[[molecule]]@bbox["x", "max"] = max(final_molecules[[molecule]]@bbox["x", "max"], 
                                                           updated_molecules[[molecule]]@bbox["x", "max"])
        final_molecules[[molecule]]@bbox["y", "max"] = max(final_molecules[[molecule]]@bbox["y", "max"], 
                                                           updated_molecules[[molecule]]@bbox["y", "max"])
      }
    }
    
    #define new offsets based on location in grid for next sample
    #go by columns, then rows
    #if modulus of next section / ncol != 1, then we are in the same row, so increase x
    #if modulus == 1, then we are in a new row, so increase y and reset x
    if((image + 1) %% n_cols == 1)
    {
      starting_x = 0
      #use boundary box so we get global maxima (don't assume similar shapes, sizes)
      starting_y = final_centroids@bbox["y", "max"] + offset 
    } else
    {
      #here we want local, not global max
      starting_x = max(updated_centroids@coords[, "x"]) + offset
    }
  }
  
  combined_fov = CreateFOV(coords = final_centroids,
                           molecules = final_molecules,
                           assay = assay,
                           key = paste0(assay, "_"))
  
  if(append == FALSE)
  {
    obj@images[1:length(obj@images)] = NULL
  }
  obj@images[[fov_name]] = combined_fov
  
  return(obj)
}


