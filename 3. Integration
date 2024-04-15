# Integrate Frontal Cortex (fc) and Motor Cortex (mc) of samples a1-a6

# Create a list of Seurat objects for fc and mc
a1_list <- list(a1_fc, a1_mc)

# Select integration features
integ_features <- SelectIntegrationFeatures(object.list = a1_list, nfeatures = 3000)

# Merge fc and mc data
a1 <- merge(a1_fc, a1_mc, merge.data = TRUE)

# Calculate percentage of mitochondrial genes
a1[["percent.mt"]] <- PercentageFeatureSet(a1, pattern = "^MT-")

# Set variable features for integration
VariableFeatures(a1) <- integ_features

# Plot violin plot
p1 <- VlnPlot(a1, c("nFeature_RNA", "nCount_RNA"), group.by = "tissue")

# Perform Harmony integration
harm <- a1
harm <- RunPCA(harm, assay = "SCT", npcs = 50, verbose = TRUE)
harm <- RunHarmony(harm, group.by.vars = "orig.ident", reduction = "pca", assay.use = "SCT", reduction.save = "harmony", verbose = TRUE)
harm <- RunUMAP(harm, reduction = "harmony", assay = "SCT", dims = 1:40, verbose = TRUE)
harm <- FindNeighbors(object = harm, reduction = "harmony", verbose = TRUE)
harm <- FindClusters(harm, resolution = 0.2, verbose = TRUE)
a1 <- harm

# Update metadata with original metadata
meta <- a1@meta.data
common_cells <- intersect(meta$cell_barcode, orig_metadata$cell_barcode)
meta$annotation_cell_class <- orig_metadata$annotation_cell_class[match(meta$cell_barcode, orig_metadata$cell_barcode)]
meta$annotation_major_cell_type <- orig_metadata$annotation_major_cell_type[match(meta$cell_barcode, orig_metadata$cell_barcode)]
meta$annotation_cell_subtype <- orig_metadata$annotation_cell_subtype[match(meta$cell_barcode, orig_metadata$cell_barcode)]
a1@meta.data <- meta

# Save Seurat object
saveRDS(a1, "seu_obj/a1.rds")
# Repeat for samples a2-a6 with fc and mc data
