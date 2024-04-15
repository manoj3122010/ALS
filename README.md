# ALS
# We analyzed sn-RNA seq data from frontal and motor cortex of six individuals suffering from c9orf72 ALS
# We compared the expression pattern of the two cortices to understand the involvement of each region and cell types

#  Some libraries required
library(Matrix)  # For working with sparse matrices
library(Seurat)  # For single-cell RNA-seq analysis
library(harmony)  # For batch correction and data integration in single-cell analysis
library(DoubletFinder)  # For doublet detection in single-cell data
library(CHETAH)  # For Cell annotation
library(MAST) #  For differential gene expression
library(ggplot2)  # For data visualization, particularly creating plots
library(dplyr) # For filtering, summarizing, and manipulating data frames
library(readr) # For reading and parsing data from CSV, TSV, and other structured data formats
library(openxlsx) # To read, write, and manipulate data in Excel format
