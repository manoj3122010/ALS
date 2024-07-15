# Define a function for gene subset based on avg_log2FC and p_val_adj
subset_genes <- function(mast_data, up_threshold, down_threshold, p_val_threshold) {
    # Subset genes based on up-regulation and p-value significance
    up_genes <- subset(mast_data, subset = avg_log2FC > up_threshold & p_val_adj < p_val_threshold)
    # Subset genes based on down-regulation and p-value significance
    down_genes <- subset(mast_data, subset = avg_log2FC < down_threshold & p_val_adj < p_val_threshold)
    # Return a list containing up-regulated and down-regulated genes
    return(list(up_genes = up_genes, down_genes = down_genes))
}

# Define thresholds for avg_log2FC based up and down-regulation, and p-value significance
up_threshold <- 0.5
down_threshold <- -0.5
p_val_threshold <- 0.05

# Subset genes for Exc_neuron, Inh_neuron, and Non_neuron for samples a1-a6
# Using lapply to apply the subset_genes function to each sample's mast_data
exc_genes <- lapply(list(a1_exc_mast, a2_exc_mast, a3_exc_mast, a4_exc_mast, a5_exc_mast, a6_exc_mast),
                     function(data) subset_genes(data, up_threshold, down_threshold, p_val_threshold))
inh_genes <- lapply(list(a1_inh_mast, a2_inh_mast, a3_inh_mast, a4_inh_mast, a5_inh_mast, a6_inh_mast),
                     function(data) subset_genes(data, up_threshold, down_threshold, p_val_threshold))
non_genes <- lapply(list(a1_non_mast, a2_non_mast, a3_non_mast, a4_non_mast, a5_non_mast, a6_non_mast),
                     function(data) subset_genes(data, up_threshold, down_threshold, p_val_threshold))

# Extract upregulated and downregulated genes for Exc_neuron, Inh_neuron, and Non_neuron for samples a1-a6
# Using lapply to extract up-regulated and down-regulated genes from each list of genes
exc_up_genes <- lapply(exc_genes, function(gene_data) gene_data$up_genes)
exc_down_genes <- lapply(exc_genes, function(gene_data) gene_data$down_genes)
inh_up_genes <- lapply(inh_genes, function(gene_data) gene_data$up_genes)
inh_down_genes <- lapply(inh_genes, function(gene_data) gene_data$down_genes)
non_up_genes <- lapply(non_genes, function(gene_data) gene_data$up_genes)
non_down_genes <- lapply(non_genes, function(gene_data) gene_data$down_genes)

# Save individual a1_exc_up files
saveRDS(exc_up_genes[[1]], "a1_exc_up.rds")
saveRDS(exc_up_genes[[2]], "a2_exc_up.rds")
saveRDS(exc_up_genes[[3]], "a3_exc_up.rds")
saveRDS(exc_up_genes[[4]], "a4_exc_up.rds")
saveRDS(exc_up_genes[[5]], "a5_exc_up.rds")
saveRDS(exc_up_genes[[6]], "a6_exc_up.rds")
# Similar saving procedure for other cell types and samples
# Example saveRDS(inh_up_genes[[4]], "a4_inh_up.rds") to save upregulated genes of sample a4 and Inh cell type

# Filter for Network Analysis
# Select upregulated genes in at least 3 samples or cell types for Exc_neuron
exc_up_genes <- c(a1_exc_up$gene, a2_exc_up$gene, a3_exc_up$gene, a4_exc_up$gene, 
  a5_exc_up$gene, a6_exc_up$gene)
gene_counts <- table(exc_up_genes)
exc_up <- names(gene_counts[gene_counts >= 3])

# Select downregulated genes in at least 3 samples or cell types for Exc_neuron
exc_down_genes <- c(a1_exc_down$gene, a2_exc_down$gene, a3_exc_down$gene, 
  a4_exc_down$gene, a5_exc_down$gene, a6_exc_down$gene)
gene_counts <- table(exc_down_genes)
exc_down <- names(gene_counts[gene_counts >= 3])

# Repeat the above process for Inh_neuron and Non_neuron
inh_up_genes <- c(a1_inh_up$gene, a2_inh_up$gene, a3_inh_up$gene, a4_inh_up$gene, 
  a5_inh_up$gene, a6_inh_up$gene)
gene_counts <- table(inh_up_genes)
inh_up <- names(gene_counts[gene_counts >= 3])

inh_down_genes <- c(a1_inh_down$gene, a2_inh_down$gene, a3_inh_down$gene, 
  a4_inh_down$gene, a5_inh_down$gene, a6_inh_down$gene)
gene_counts <- table(inh_down_genes)
inh_down <- names(gene_counts[gene_counts >= 3])

non_up_genes <- c(a1_non_up$gene, a2_non_up$gene, a3_non_up$gene, a4_non_up$gene, 
  a5_non_up$gene, a6_non_up$gene)
gene_counts <- table(non_up_genes)
non_up <- names(gene_counts[gene_counts >= 3])

non_down_genes <- c(a1_non_down$gene, a2_non_down$gene, a3_non_down$gene, 
  a4_non_down$gene, a5_non_down$gene, a6_non_down$gene)
gene_counts <- table(non_down_genes)
non_down <- names(gene_counts[gene_counts >= 3])

# Save the lists of upregulated and downregulated genes for further analysis
saveRDS(exc_up, "exc_up_genes.rds")
saveRDS(exc_down, "exc_down_genes.rds")
saveRDS(inh_up, "inh_up_genes.rds")
saveRDS(inh_down, "inh_down_genes.rds")
saveRDS(non_up, "non_up_genes.rds")
saveRDS(non_down, "non_down_genes.rds")
# Filter for ALS genes
# `als_genes` is a tibble containing ALS-related genes obtained from the OpenTargets platform.
# It includes attributes such as the gene symbol and overall association score with ALS.
# Additionally, it provides information on genetic associations, somatic mutations, drugs targeting the genes,
# pathways/systems biology relevance, text mining evidence, RNA expression, and more.

# Combine upregulated genes from Exc_neuron, Inh_neuron, and Non_neuron into a single list
up <- c(exc_up_genes, inh_up_genes, non_up_genes)

# Create a data frame with the upregulated genes
up <- data.frame(symbol = up)

# Filter ALS genes that are upregulated
als_up <- als_genes %>%
     semi_join(up, by = "symbol")

# Combine downregulated genes from Exc_neuron, Inh_neuron, and Non_neuron into a single list
down <- c(exc_down_genes, inh_down_genes, non_down_genes)

# Create a data frame with the downregulated genes
down <- data.frame(symbol = down)

# Filter ALS genes that are downregulated
als_down <- als_genes %>%
     semi_join(down, by = "symbol")

# Find common genes between upregulated and downregulated ALS genes
als_common <- intersect(als_up$symbol, als_down$symbol)

# Create a data frame with the common genes
common <- data.frame(symbol = als_common)

# Filter ALS genes that are common between upregulated and downregulated
als_common <- als_genes %>%
     semi_join(common, by = "symbol")

# Write the filtered ALS genes to Excel files
write.xlsx(als_up, file = "als_up.xlsx", row.names = FALSE)
write.xlsx(als_down, file = "als_down.xlsx", row.names = FALSE)
write.xlsx(als_common, file = "als_common.xlsx", row.names = FALSE)
