# Doublet Detection for Frontal Cortex (fc) and Motor Cortex (mc) in samples a1-a6

# Process fc data for sample a1
a1_fc <- SCTransform(a1_fc)  # Perform SCTransform
a1_fc <- RunPCA(a1_fc)  # Run PCA
a1_fc <- RunUMAP(a1_fc, dims = 1:30)  # Run UMAP
ElbowPlot(a1_fc)  # Plot Elbow Plot
sweep.res.list_a1_fc <- paramSweep(a1_fc, PCs = 1:20, sct = TRUE)  # Parameter sweep
sweep.stats_a1_fc <- summarizeSweep(sweep.res.list_a1_fc, GT = FALSE)  # Summarize sweep results
bcmvn_a1_fc <- find.pK(sweep.stats_a1_fc)  # Find optimal pK
nExp_poi_a1_fc <- round(0.053 * nrow(a1_fc@meta.data))  # Calculate expected number of doublets
a1_fc_plot <- ggplot(bcmvn_a1_fc, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()  # Create plot
ggsave("a1_fc_pK.png", a1_fc_plot, height = 4, width = 10, dpi = 750)  # Save plot
a1_fc <- doubletFinder(a1_fc, PCs = 1:20, pN = 0.25, pK = 0.19, nExp = nExp_poi_a1_fc, reuse.pANN = FALSE, sct = TRUE)  # Find doublets

# Process mc data for sample a1
a1_mc <- SCTransform(a1_mc)  # Perform SCTransform
a1_mc <- RunPCA(a1_mc)  # Run PCA
a1_mc <- RunUMAP(a1_mc, dims = 1:30)  # Run UMAP
ElbowPlot(a1_mc)  # Plot Elbow Plot
sweep.res.list_a1_mc <- paramSweep(a1_mc, PCs = 1:20, sct = TRUE)  # Parameter sweep
sweep.stats_a1_mc <- summarizeSweep(sweep.res.list_a1_mc, GT = FALSE)  # Summarize sweep results
bcmvn_a1_mc <- find.pK(sweep.stats_a1_mc)  # Find optimal pK
nExp_poi_a1_mc <- round(0.053 * nrow(a1_mc@meta.data))  # Calculate expected number of doublets
a1_mc_plot <- ggplot(bcmvn_a1_mc, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()  # Create plot
ggsave("a1_mc_pK.png", a1_mc_plot, height = 4, width = 10, dpi = 750)  # Save plot
a1_mc <- doubletFinder(a1_mc, PCs = 1:20, pN = 0.25, pK = 0.08, nExp = nExp_poi_a1_mc, reuse.pANN = FALSE, sct = TRUE)  # Find doublets

# Repeat for samples a2-a6 with fc and mc data

# Add metadata for sample a1
meta_a1_fc <- a1_fc@meta.data  # Get metadata for a1_fc
colnames(meta_a1_fc)[colnames(meta_a1_fc) == "DF.classifications_0.25_0.19_854"] <- "Doublet"  # Rename classification column
colnames(meta_a1_fc)[colnames(meta_a1_fc) == "pANN_0.25_0.19_854"] <- "pANN_dbfinder"  # Rename pANN column
meta_a1_fc$cell_barcode <- rownames(meta_a1_fc)  # Add cell barcode column
meta_a1_fc$sample <- "a1"  # Add sample column
meta_a1_fc$tissue <- "fc"  # Add tissue column
meta_a1_mc <- a1_mc@meta.data  # Get metadata for a1_mc
colnames(meta_a1_mc)[colnames(meta_a1_mc) == "DF.classifications_0.25_0.08_472"] <- "Doublet"  # Rename classification column
colnames(meta_a1_mc)[colnames(meta_a1_mc) == "pANN_0.25_0.08_472"] <- "pANN_dbfinder"  # Rename pANN column
meta_a1_mc$cell_barcode <- rownames(meta_a1_mc)  # Add cell barcode column
meta_a1_mc$sample <- "a1"  # Add sample column
meta_a1_mc$tissue <- "mc"  # Add tissue column
a1_fc@meta.data <- meta_a1_fc  # Update metadata for a1_fc
a1_mc@meta.data <- meta_a1_mc  # Update metadata for a1_mc
saveRDS(a1_fc, "seu_obj/a1_fc.rds")  # Save Seurat object for a1_fc
saveRDS(a1_mc, "seu_obj/a1_mc.rds")  # Save Seurat object for a1_mc

# Repeat for samples a2-a6 with fc and mc data
