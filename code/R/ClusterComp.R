library(Seurat)
library(patchwork)

ClusterComp <- function (seurat_object){
  # Get the cluster identities and original identities
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
    labs(x = "Cluster", y = "Cluster Composition (Cell Count)", fill = "Sample", title = "Sample Composition of Plasma clusters")
  
  fraction_plot <- ggplot(df_counts, aes(x = Cluster, y = Fraction, fill = OrigIdent)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(x = "Cluster", y = "Cluster Composition (%)", fill = "Sample", title = "Sample Composition of Plasma clusters")
  
  count_plot + fraction_plot
}