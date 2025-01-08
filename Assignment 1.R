#Q1-1
library(monocle3)
cds <- readRDS("~/Documents/Columbia/Sem 3/Funct Genomics/Dataset1/dataset1_final.rds")
expr_matrix <- exprs(cds)

#Dimensions
dimensions <- dim(expr_matrix)
cat("The expression matrix has", dimensions[1], "genes and", dimensions[2], "cells.\n")

#Q1-b
library(ggplot2)
library(readr)

#Number of cells
cell_meta_df <- read_csv('~/Documents/Columbia/Sem 3/Funct Genomics/Dataset1/dataset1_final_cellmeta.csv')
num_cells <- nrow(cell_meta_df)

#Median number of UMIs per cell
median_umis <- median(cell_meta_df$n.umi)
print(paste("Total number of cells:", num_cells))
print(paste("Median number of UMIs per cell:", median_umis))

#Violin plot with boxplot overlay and median line
ggplot(cell_meta_df, aes(x = 1, y = n.umi)) +  # Set x to 1 as a dummy variable
  geom_violin(fill = "lightblue") +
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA) +
  geom_hline(yintercept = median_umis, color = "red", linetype = "dashed") +
  labs(title = "Distribution of UMIs per Cell", y = "Number of UMIs", x = "") +
  theme_minimal() +
  annotate("text", x = 1.2, y = median_umis + 200, label = paste("Median:", round(median_umis, 1)), color = "red")

#Q1-c
cds <- readRDS('~/Documents/Columbia/Sem 3/Funct Genomics/Dataset1/dataset1_final.rds')
cds <- detect_genes(cds)
num_genes <- nrow(rowData(cds))

#Number of genes
genes_per_cell <- colData(cds)$num_genes_expressed

#Median number of genes quantified per cell
median_genes_per_cell <- median(genes_per_cell)
print(paste("Total number of genes:", num_genes))
print(paste("Median number of genes per cell:", median_genes_per_cell))

#Violin plot with boxplot overlay and median line for genes per cell
ggplot(data = as.data.frame(genes_per_cell), aes(x = 1, y = genes_per_cell)) +
  geom_violin(fill = "lightblue") +
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA) +
  geom_hline(yintercept = median_genes_per_cell, color = "red", linetype = "dashed") +
  labs(title = "Distribution of Quantified Genes per Cell", y = "Number of Genes Quantified", x = "") +
  theme_minimal() +
  annotate("text", x = 1.2, y = median_genes_per_cell + 200, label = paste("Median:", round(median_genes_per_cell, 1)), color = "red")

#Q1-d
#Column names
colnames(colData(cds))
unique_perturbations <- unique(colData(cds)$crispr_target)

#Total number of cells for each genetic perturbation
summary_table <- table(colData(cds)$crispr_target)
summary_df <- as.data.frame(summary_table)
colnames(summary_df) <- c("Genetic_Perturbation", "Number_of_Cells")
print(summary_df)

#Q1-d-b
#Genetic perturbations
valid_perturbations <- c("CHEK1", "HNRNPC", "non-targeting", "RUVBL1", "SSRP1", "SUPT5H", "SUPT6H")

#Filter dataset to valid perturbations
filtered_data <- colData(cds)[colData(cds)$crispr_target %in% valid_perturbations, ]

#Violin plot for distribution of UMIs per cell by genetic perturbation
ggplot(as.data.frame(filtered_data), aes(x = crispr_target, y = n.umi)) +
  geom_violin(fill = "lightblue") +
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", shape = 20, size = 3, color = "red", show.legend = FALSE) +
  labs(title = "Distribution of UMIs per Cell by Genetic Perturbation", y = "Number of UMIs", x = "Genetic Perturbation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Q1-d-c
#Number of genes (rows in rowData) and median genes per cell (colData)
num_genes <- nrow(rowData(cds))
genes_per_cell <- colData(cds)$num_genes_expressed
median_genes_per_cell <- median(genes_per_cell)

cat("Total number of genes:", num_genes, "\n")
cat("Median number of genes per cell:", median_genes_per_cell, "\n")

#Violin plot for number of genes per cell by genetic perturbation
library(ggplot2)

ggplot(as.data.frame(colData(cds)), aes(x = crispr_target, y = num_genes_expressed)) +
  geom_violin(fill = "lightblue") +
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA) +
  geom_hline(yintercept = median_genes_per_cell, color = "red", linetype = "dashed") +
  labs(title = "Distribution of Genes Quantified per Cell by Genetic Perturbation", 
       y = "Number of Genes Quantified", x = "Genetic Perturbation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Q1-d