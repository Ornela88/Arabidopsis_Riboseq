# This script calculates the correlation between replicates of Time point 1 after the calculation of the relative coverage. 

rm(list = ls())


# PRELIMINARIES =======================================

#sets the working directory
#setwd("~/rds/rds-dc702-bio2_core_rds/Personal_folders_int/dc702/Cooperative_work_Ornela_Sophia/Davide_Ribo_mtx")

#uploads the required packages
require(data.table)

# DATA UPLOAD =======================================
# Uploads the matrices coverages

matrix_cov_A1 <-fread("./matrix_cov_A1_1.csv", nThread = 20, stringsAsFactors = F)

matrix_cov_A1<- as.data.frame(matrix_cov_A1)
rownames(matrix_cov_A1) <- matrix_cov_A1[,1]
matrix_cov_A1<-matrix_cov_A1[,-1]
matrix_cov_A1<-as.matrix(matrix_cov_A1)


matrix_cov_A2 <-fread("./matrix_cov_A1_2.csv", nThread = 20, stringsAsFactors = F)

matrix_cov_A2<- as.data.frame(matrix_cov_A2)
rownames(matrix_cov_A2) <- matrix_cov_A2[,1]
matrix_cov_A2<-matrix_cov_A2[,-1]
matrix_cov_A2<-as.matrix(matrix_cov_A2)

matrix_cov_A3 <-fread("./matrix_cov_A1_3.csv", nThread = 20, stringsAsFactors = F)

matrix_cov_A3<- as.data.frame(matrix_cov_A3)
rownames(matrix_cov_A3) <- matrix_cov_A3[,1]
matrix_cov_A3<-matrix_cov_A3[,-1]
matrix_cov_A3<-as.matrix(matrix_cov_A3)

# identifies the transcripts in common among the replicates
transcripts_in_common <- Reduce(intersect, list(rownames(matrix_cov_A1),
                                                rownames(matrix_cov_A2),
                                                rownames(matrix_cov_A3)))


# produces matrices coverage with only the transcripts in common
matrix_cov_A1_common <- matrix_cov_A1[which(rownames(matrix_cov_A1) %in%  transcripts_in_common ), ]
matrix_cov_A2_common <- matrix_cov_A2[which(rownames(matrix_cov_A2) %in%  transcripts_in_common ), ]
matrix_cov_A3_common <- matrix_cov_A3[which(rownames(matrix_cov_A3) %in%  transcripts_in_common ), ]

# computes the relative coverage
compute_relative_coverage <- function(x) (x/sum(x, na.rm = T))

matrix_cov_A1_common_relative <- t(apply(matrix_cov_A1_common, 1, compute_relative_coverage))
matrix_cov_A2_common_relative <- t(apply(matrix_cov_A2_common, 1, compute_relative_coverage))
matrix_cov_A3_common_relative <- t(apply(matrix_cov_A3_common, 1, compute_relative_coverage))

# computes the lengths of the common transcripts
compute_transcripts_lengths <- function(x) length(which(!is.na(x)))

transcripts_lengths_A1 <- apply( matrix_cov_A1_common, 1, compute_transcripts_lengths )
transcripts_lengths_A2 <- apply( matrix_cov_A2_common, 1, compute_transcripts_lengths )
transcripts_lengths_A3 <- apply( matrix_cov_A3_common, 1, compute_transcripts_lengths )


# select from the coverage matrix the raw that you want to calculate for the correlation

A1_2 <- matrix_cov_A2_common_relative["AT1G04290.2",c(1:transcripts_lengths_A1["AT1G04290.2"])]
A1_1 <- matrix_cov_A1_common_relative["AT1G04290.2",c(1:transcripts_lengths_A1["AT1G04290.2"])]
A1_3 <- matrix_cov_A3_common_relative["AT1G04290.2",c(1:transcripts_lengths_A1["AT1G04290.2"])]

# produces matrices coverage considering 

M<- cbind(A1_1,A1_2,A1_3) 


# CORRELATION ANALYSIS AND SCATTERPLOTS =======================================

#Pearson's Correlation for Time Point 1 replicate 1 vs 2
cor_1_2 <-cor.test(M[,1],M[,2]) 
cor_1_2 
# Creating the plot
x <- M[,1]
y <- M[,2]
pdf("A1_cor_1vs2_AT1G04290.2.pdf")
plot(x, y, main = "Pearson's Correlation TP1 1vs2 for AT1G04290.2 ",
     xlab = "TP1 rep 1 ", ylab = "TP1 rep 2",
     pch = 16, sub ="r = 0.078")
abline(lm(y ~ x), col = "blue")
dev.off()

#Pearson's Correlation for Time Point 1 replicate 1 vs 3
cor_1_3 <-cor.test(M[,1],M[,3])
cor_1_3 
# Creating the plot
x <- M[,1]
y <- M[,3]
pdf("A1_cor_1vs3_AT1G04290.2.pdf")
plot(x, y, main = "Pearson's Correlation TP1 1vs3 AT1G04290.2 ",
     xlab = "TP1 rep 1", ylab = "TP1 rep 3",
   pch = 16, sub = "r = -0.14")
abline(lm(y ~ x), col = "blue")
dev.off()

#Pearson's Correlation for Time Point 1 replicate 1 vs 3

cor_2_3 <-cor.test(M[,2],M[,3]) 
cor_2_3
# Creating the plot
x <-M[,2]
y <-M[,3]
pdf("A1_cor_2vs3_AT1G04290.2.pdf")
plot(x, y, main = "Pearson's Correlation TP1 2vs3 for AT1G04290.2",
     xlab = "TP1 rep 2 ", ylab = "TP1 rep 3",
  pch = 16, sub = "r = 0.21 ")
abline(lm(y ~ x), col = "blue")
dev.off()

