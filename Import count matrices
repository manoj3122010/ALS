# Read spliced data for samples a1-a6 frontal cortex (fc) and motor cortex (mc)
a1_fc <- ReadMtx(mtx = "a1/fc/spliced.mtx", 
                 features = "a1/fc/genes.txt", 
                 cells = "a1/fc/spliced.barcodes.txt", 
                 feature.column = 1, 
                 mtx.transpose = TRUE)
a1_mc <- ReadMtx(mtx = "a1/mc/spliced.mtx", 
                 features = "a1/mc/genes.txt", 
                 cells = "a1/mc/spliced.barcodes.txt", 
                 feature.column = 1, 
                 mtx.transpose = TRUE)

# Repeat for a2, a3, a4, a5, a6 with both fc and mc reads

# Create Seurat objects for each sample's fc and mc data
a1_fc <- CreateSeuratObject(a1_fc, project = "a1_fc", min.cells = 5, min.features = 500)
a1_mc <- CreateSeuratObject(a1_mc, project = "a1_mc", min.cells = 5, min.features = 500)

# Repeat for a2, a3, a4, a5, a6 with both fc and mc reads

# Subset cells based on RNA feature counts for each sample's fc and mc data
a1_fc <- subset(a1_fc, subset = nFeature_RNA >= 500 & nFeature_RNA <= 8000)
a1_mc <- subset(a1_mc, subset = nFeature_RNA >= 500 & nFeature_RNA <= 8000)

# Repeat for a2, a3, a4, a5, a6 with both fc and mc reads

# Subset cells based on RNA molecule counts for each sample's fc and mc data
a1_fc <- subset(a1_fc, subset = nCount_RNA >= 700 & nCount_RNA <= 20000)
a1_mc <- subset(a1_mc, subset = nCount_RNA >= 700 & nCount_RNA <= 20000)

# Repeat for a2, a3, a4, a5, a6 with both fc and mc reads

# Generate Violin Plots for each sample's fc and mc data
VlnPlot(a1_fc, c("nFeature_RNA", "nCount_RNA"))  # Violin Plot for a1_fc
VlnPlot(a1_mc, c("nFeature_RNA", "nCount_RNA"))  # Violin Plot for a1_mc

# Repeat for a2, a3, a4, a5, a6 with both fc and mc reads
