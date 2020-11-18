#library
options(stringsAsFactors = F)
#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("TxDb.Athaliana.BioMart.plantsmart28")
#BiocManager::install("Rsamtools")
library(Rsamtools)
library(RiboProfiling)
library(TxDb.Athaliana.BioMart.plantsmart28)
library(org.At.tair.db)
library("GO.db")
library("GOstats")
library(GenomicAlignments)
library(xlsx)
library(BSgenome)
library(BSgenome.Athaliana.TAIR.TAIR9)
library(Biostrings)
library(GenomeInfoDb)
library(BSgenome)
library('GenomicFeatures')
library(GenomeInfoDb)
library(GenomicRanges)
library(Biostrings)
library(GenomicFeatures)
library(ggplot2)
library(grid)
library(systemPipeR)
library(S4Vectors)
library(GGally)
library(RiboseQC)
library(limma)
library(Rsubread)
library(edgeR)
library(biomaRt)
library(RColorBrewer)
library(DESeq2, quietly = TRUE)
library(ape, warn.conflicts = FALSE)
library(gplots)

######Read multiple BAM files together
WORKDIR="/usr/data/bgfs1/maloku/riboseq_data_Ornela/bam_files/";
dir.create(WORKDIR)
setwd(WORKDIR)
WORK_RiboProfiling ="/usr/home/maloku/RiboProfiling_test/"; 
dir.create(WORK_RiboProfiling_Result)
BAI <- list.files(WORKDIR, pattern = "bam.bai$",full.names = T)
SAMPLES_bai <- gsub(".bai","",list.files("/usr/data/bgfs1/maloku/riboseq_data_Ornela/bam_files/", pattern = ".bai$"))
WORK_RiboProfiling_Result <-paste0(WORK_RiboProfiling, "Results/")
BAMS <- list.files(WORKDIR, pattern = ".bam$",full.names = T)

################### FeatureCounts to quantifying read counts for each gene
fc <-featureCounts(files=BAMS, annot.ext="/usr/home/maloku/RiboProfiling_test/Arabidopsis_thaliana.TAIR10.45.gtf",
                   useMetaFeatures = F, isPairedEnd = F, isGTFAnnotationFile = T, nthreads = 20, annot.inbuilt="TAIR10")
names(fc)

###Normalization step

# Convert counts to dge object

d0 <- DGEList(fc$counts) 
d0$genes <- fc$annotation[,c( "GeneID","Length"), drop=FALSE]
d0 <-calcNormFactors(d0)
##### or
isexpr <- rowSums(cpm(d0) > 2) >= 3 ## Create filter genes that have less than 2 counts per million in at least three libraries
d0<- d0[isexpr,]
d0.rpkm <- rpkm(d0,d0$genes$Length) # create object with rpkm values estimated from the counts
write.csv(file = '/usr/home/maloku/RiboProfiling_test/CountValues.csv', d0$counts) # Save raw count data
write.csv(file = '/usr/home/maloku/RiboProfiling_test/RPKMValues.csv', d0.rpkm) # Save RPKM values


## definition of a model matrix. This design matrix simply links each group to the samples that belong to it.
targets <-readTargets("/usr/home/maloku/RiboProfiling_test/samples.csv",sep=",")
time <-targets$time
mm <- model.matrix(~0+time)
colnames(time) <- levels(time)
mm



##############Differential expression using voom ###################
# use this when the library sizes are quite variable between samples, then the voom approach is theoretically more powerful than limma-trend

y <- voom(d0, mm, plot=T, normalize ="quantile")
names(y)
cor(y$E)

#correlations among samples and graphically show similarity among samples with a MDS plot.
write.csv(cor(y$E), file="/usr/home/maloku/RiboProfiling_test/d0_correlation_yE.csv")
plotMDS(y)
par(mfrow=c(1,2))
boxplot(y$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
#blue horizontal line that corresponds to the median logCPM
abline(h=median(y$E),col="blue")

#### fitting  linear models

fitV <-lmFit(y,mm)
names(fitV)
head(coef(fitV))
write.csv(fitV, file="/usr/home/maloku/RiboProfiling_test/fitV_model_all.csv")
fitV <- eBayes(fitV)
#Summarize the number of significant genes, p < 0.05 after BH-correction
top_table_all <-topTable(fitV, coef=ncol(mm), adjust.method="BH",sort.by="none")


#TestResults matrix. This is numeric matrix of 0's, 1's and -1's indicating significance of a test or membership of a set.
results <- decideTests(fitV)
#write.csv(vennCounts(results, include="both"), file="/usr/home/maloku/RiboProfiling_test/fitV_model_up_and_down_genes.csv")

# Venn diagram that shows the sets of sign. genes that are common among the contrasts
vennCounts(results)
vennDiagram(results[,4:7], circle.col=c("red","green","blue","yellow"), main = 'All significant below adjusted P-value 0.05')
vennDiagram(results[,4:7], include = 'down', circle.col=c("red","green","blue","yellow"), main = 'Down-regulated below adjusted P-value 0.05')
vennDiagram(results[,4:7], include = 'up', circle.col=c("red","green","blue","yellow"), main = 'Up-regulated below adjusted P-value 0.05')

##comparison between groups and contrast estimation
contr <-makeContrasts(T1vsT3=timeT1-timeT3,T5vsT7=timeT5-timeT7,T9vsT11=timeT9-timeT11,T13vsT15=timeT13-timeT15,
                      T17vsT19=timeT17-timeT19,T21vsT23=timeT21-timeT23,
                      levels=colnames(coef(fitV)))
tmp <- contrasts.fit(fitV,contr)
tmp <-eBayes(tmp)
summa.fit <- decideTests(tmp)
summary(summa.fit)
top.table <-topTable(tmp,sort.by="none", adjust.method="BH")
head(top.table,5)
write.table(top.table, file ="/usr/home/maloku/RiboProfiling_test/top_table_limma.txt", 
                   row.names=F, sep="\t", quote=F)

#### Plot P value and  adj P value
pdf(paste0(WORK_RiboProfiling_Result, "P_VAl_adj_P_val.pdf"), width = 15, height = 15)
layout(matrix(1:2, ncol=1))
tab1<-cumsum(table(cut(top.table$P.Value, breaks=c(0,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.05, 1e-1,1))))
barplot(tab1, cex.names=0.7, col="darkolivegreen2", space=0, main="Cumulative number of genes attaining defined significance - unadjusted P values Correlation Analyses", xlab="Gene significance", ylab="Gene number")
for(i in 1:9)text(x=i-0.5, y=1, tab1[i], pos=3)
tab4<-cumsum(table(cut(top.table$adj.P.Val, breaks=c(0,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.05, 1e-1,1))))
barplot(tab1, cex.names=0.7, col="darkolivegreen2", space=0, main="Cumulative number of genes attaining defined significance - adjusted P values Correlation Analyses", xlab="Gene significance", ylab="Gene number")
for(i in 1:9)text(x=i-0.5, y=1, tab4[i], pos=3)
dev.off()

################comparison between groups and contrast estimation for T1 vs T3
contrT1_T3 <-makeContrasts(T1vsT3=timeT1-timeT3,
                      levels=colnames(coef(fitV)))
tmp1 <- contrasts.fit(fitV,contrT1_T3)
tmp1 <-eBayes(tmp1)
top.table1 <-topTable(tmp1,sort.by="none",adjust.method="BH")
write.csv(top.table1, file ="/usr/home/maloku/RiboProfiling_test/topT1_T3_limma.csv")

#### Plot P value and  adj P value for T1 vs T3

pdf(paste0(WORK_RiboProfiling_Result, "P_VAl_adj_P_val_T1_T3.pdf"), width = 15, height = 15)
layout(matrix(1:2, ncol=1))
tab1<-cumsum(table(cut(top.table1$P.Value, breaks=c(0,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.05, 1e-1,1))))
barplot(tab1, cex.names=0.7, col="darkolivegreen2", space=0, main="Cumulative number of genes attaining defined significance - unadjusted P values Correlation Analyses", xlab="Gene significance", ylab="Gene number")
for(i in 1:9)text(x=i-0.5, y=1, tab1[i], pos=3)
tab4<-cumsum(table(cut(top.table1$adj.P.Val, breaks=c(0,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.05, 1e-1,1))))
barplot(tab1, cex.names=0.7, col="darkolivegreen2", space=0, main="Cumulative number of genes attaining defined significance - adjusted P values Correlation Analyses", xlab="Gene significance", ylab="Gene number")
for(i in 1:9)text(x=i-0.5, y=1, tab4[i], pos=3)
dev.off()

######scatterplot of T1 effect vs T3 effect

T1<- topTable(fitV[,1], number = Inf, sort.by = "t", coef = grep("timeT1", 
                                                                 colnames(coef(fitV[,1]))))
T3 <- topTable(fitV, number = Inf, sort.by = "t", coef = grep("timeT3", 
                                                                         colnames(coef(fitV))))
names(T1)
names(T3)
smoothScatter(T1$t ~ T3$t, xlim = c(-20, 20), ylim = c(-20, 20), xlab = "t-statistic for effect of T1", 
              ylab = "t-statistic for effect of T3")
abline(a = 0, b = 1, col = "orange")


##################heatmap plot
png("/usr/home/maloku/RiboProfiling_test/heatmaps_in_r",  width=700, height=700)
logCPM <- cpm(d0, log=TRUE, prior.count=3)
var_genes <- apply(logCPM, 1, var)
select_var <- names(sort(var_genes, decreasing=TRUE))[1:20]
my_palette <- colorRampPalette(c("red", "yellow", "green")) 
highly_variable_lcpm <- logCPM[select_var,]
heatmap.2(highly_variable_lcpm,col=rev(my_palette(70)),trace="none", main="Top 300 most variable genes across samples",scale="row")

###or use pheatmap 
pheatmap(highly_variable_lcpm, scale = "row", clustering_distance_rows = "correlation", 
         clustering_distance_cols = "correlation")
dev.off()



