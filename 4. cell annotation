# Cell Annotation using CHETAH Classifier

# For reference the expression data from original publication was downloaded and used
# Get reference cell annotations and counts
ref_ct <- ref$annotation_cell_class
ref_counts <- GetAssayData(ref, layer = 'counts')

# Get input PCA and counts
input_pca <- a1@reductions$pca
input_counts <- GetAssayData(a1, layer = 'counts')

# Create SingleCellExperiment object for reference
reference <- SingleCellExperiment(assays = list(counts = ref_counts), colData = DataFrame(celltypes = ref_ct))

# Create SingleCellExperiment object for input with reduced dimensions
input <- SingleCellExperiment(assays = list(counts = input_counts), reducedDims = SimpleList(TSNE = input_pca))

# Run CHETAH classifier
a1_cell_chetah <- CHETAHclassifier(input = input, ref_cells = reference)

# Check the distribution of cell types predicted by CHETAH
table(as.data.frame(input$celltype_CHETAH))

# Identify common barcodes between input and CHETAH results
common_barcodes <- intersect(colnames(a1), colnames(a1_cell_chetah))

# Assign CHETAH cell types to the input Seurat object
a1$celltype_CHETAH <- a1_cell_chetah$celltype_CHETAH[common_barcodes]

#repeat for samples a2 - a6
