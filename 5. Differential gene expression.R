# Subset cells based on CHETAH cell types for sample a1
a1_exc <- subset(a1, subset = celltype_CHETAH == "Exc_neuron")
a1_inh <- subset(a1, subset = celltype_CHETAH == "Inh_neuron")
a1_non <- subset(a1, subset = celltype_CHETAH == "Non_neuron")

# Assign identities based on tissue type
Idents(a1_exc) <- a1_exc$tissue
Idents(a1_inh) <- a1_inh$tissue
Idents(a1_non) <- a1_non$tissue

# Prepare for SCTFindMarkers
a1_exc <- PrepSCTFindMarkers(a1_exc)
a1_inh <- PrepSCTFindMarkers(a1_inh)
a1_non <- PrepSCTFindMarkers(a1_non)

# Perform Wilcoxon rank sum test for differential expression
a1_exc_mast <- FindMarkers(a1_exc, test.use = "MAST", ident.1 = "fc", ident.2 = "mc")
a1_inh_mast <- FindMarkers(a1_inh, test.use = "MAST", ident.1 = "fc", ident.2 = "mc")
a1_non_mast <- FindMarkers(a1_non, test.use = "MAST", ident.1 = "fc", ident.2 = "mc")

# Add sample and cell type information
a1_exc_mast$sample <- "a1"
a1_inh_mast$sample <- "a1"
a1_non_mast$sample <- "a1"
a1_exc_mast$celltype <- "Exc_neuron"
a1_inh_mast$celltype <- "Inh_neuron"
a1_non_mast$celltype <- "Non_neuron"
a1_exc_mast$gene <- rownames(a1_exc_mast)
a1_inh_mast$gene <- rownames(a1_inh_mast)
a1_non_mast$gene <- rownames(a1_non_mast)

# Save differential expression results for sample a1
saveRDS(a1_exc_mast, "a1/deg/a1_exc_mast.rds")
saveRDS(a1_inh_mast, "a1/deg/a1_inh_mast.rds")
saveRDS(a1_non_mast, "a1/deg/a1_non_mast.rds")

# Repeat the above process for samples a2-a6
