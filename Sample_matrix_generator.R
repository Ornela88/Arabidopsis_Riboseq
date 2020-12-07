# PRELIMINARIES ============================================================

rm(list = ls())

# Sets the work directory
#setwd("~/Google Drive File Stream/My Drive/Ribosome_reproducibility_SAILab/Progetto_Giorgia/Reproducibility_E.coli/Comparisons_Python_vs_R/Small_example_for_debugging/")

# UPLOADS AND MANIPULATES THE INPUT FILES ============================================================

# Uploads the small test bed file 
bed_table_tair10 <- as.data.frame(read.table("/usr/home/maloku/Desktop/Davide_Ribo/RFP_A1_1B.bed", sep="\t"))

# Picks the relevant columns from the bed_table_FP
bed_table_tair10_reduced <- bed_table_tair10[,c(1,2,3)]

# Produces a vector containing the lengths of the reads, computed from the bed_table_FP_reduced

bed_table_tair10_reduced$reads_lengths <- bed_table_tair10_reduced[,3] - bed_table_tair10_reduced[,2] + 1

# Creates a list of the genes contained in bed_table_FP_reduced

genes_bed_tair10 <- unique(bed_table_tair10_reduced$V1)

# Uploads a file containing the list of all the genes of Athaliana and their length

all_genes_lengths_tair10 <- as.data.frame(read.csv("/usr/home/maloku/Desktop/Davide_Ribo/tair10_genes_lengths.csv", sep=","),stringsAsFactors=FALSE)

# Consider only the entries of the all_genes_lengths table that have a correspondent entry in the genes_bed vector

genes_bed_lengths_tair10 <- all_genes_lengths_tair10[which(all_genes_lengths_tair10$genenameonly %in% genes_bed_tair10),]
genes_bed_lengths_tair10 <- genes_bed_lengths_tair10[,-1]


# PRODUCES THE MATRIX COVERAGE =============================================================
# The matrix coverage has the ORFs (or genes) names in teh row and the positions ov the single nucleotides in the columns

# Initialises the matrix
matrix_coverage <- matrix(nrow=nrow(genes_bed_lengths_tair10), ncol=max(all_genes_lengths_tair10$genelenghthonly))
rownames(matrix_coverage) <- genes_bed_lengths_tair10$genenameonly

# Produces the matrix, row by row 
for (i in 1 : nrow(matrix_coverage) ){
  
  sub_bed_table <- bed_table_tair10_reduced[which(bed_table_tair10_reduced$V1 %in% genes_bed_lengths_tair10$genenameonly[i]),]
  
  vector_coverage <- vector(mode="numeric", length=genes_bed_lengths_tair10$genelenghthonly[i])
  
  for(l in 1:nrow(sub_bed_table)){
    
    vector_coverage[seq(sub_bed_table[l,2]+1,sub_bed_table[l,3]+1)] <- vector_coverage[seq(sub_bed_table[l,2]+1,sub_bed_table[l,3]+1)] + 1
  }
  
  matrix_coverage[i,1:length(vector_coverage)] <- vector_coverage    
  rownames(matrix_coverage)[i] <- as.vector(genes_bed_lengths_tair10$genenameonly)[i]
}

relative_counts <- function(x) x/sum(x, na.rm = T)
matrix_coverage_relative <- t(apply(matrix_coverage, 1, relative_counts))


