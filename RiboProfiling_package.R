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
#######RiboProfiling#####
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
library(reshape)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)
######Note that TAIR10 is an "annotation release" based on the same genome assembly as TAIR9####
######https://davetang.org/muse/2017/08/08/getting-started-arabidopsis-thaliana-genomics/
txdb <- TxDb.Athaliana.BioMart.plantsmart28
#columns(txdb)
#all_gene <- genes(txdb)
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

####Upload bam files with BamFile ###############################################

A1_3_bam <-BamFile("/usr/data/bgfs1/maloku/riboseq_data_Ornela/bam_files/RFP_A1_3B_genome_Aligned.sortedByCoord.out.bam")
A3_3_bam <-BamFile("/usr/data/bgfs1/maloku/riboseq_data_Ornela/bam_files/RFP_A3_3B_genome_Aligned.sortedByCoord.out.bam")

######load the BAM file and index files

listInputBam <-c(A1_3_bam,A3_3_bam)
png(filename="/usr/home/maloku/RiboProfiling_test/bam_results/covData.Rd_%03d_prova.png",width=550, height=660)
covData <- riboSeqFromBAM(listInputBam,
                        txdb=txdb,
                        param=ScanBamParam(),
                        listShiftValue())
                       
dev.off()

#get all CDSs by transcript
cds <- GenomicFeatures::cdsBy(txdb, by="tx", use.names=TRUE)
#get all exons by transcript
exonGRanges <- GenomicFeatures::exonsBy(txdb, by="tx", use.names=TRUE)
#get the per transcript relative position of start and end codons
cdsPosTransc <- orfRelativePos(cds, exonGRanges)

############ Create a GAlignments from a BAM using readGAlignment
aln <- readGAlignments(A3_3_bam,param = NULL,
                       use.names=FALSE,
                       with.which_label =FALSE )

################################plot the distribution of match length

png(filename="/usr/home/maloku/RiboProfiling_test/Results/matchLenDistr.Rd_%03d_a1.png",width=550, height=660)
matchLenDistr<- histMatchLength(aln, 1)
matchLenDistr[[2]]
dev.off()


########add start coverage plot around the TSS and determin the offset based on the start of the read

alnGRanges <- readsToStartOrEnd(aln, what="start")
oneBinRanges <- aroundPromoter(txdb, alnGRanges, percBestExpressed=0.03,flankSize = 20)
listPromoterCov <-readStartCov(alnGRanges,oneBinRanges,
                               matchSize=c(27:32),
                               fixedInterval=c(-20,20), renameChr="aroundTSS", charPerc="perc")

#View(as.data.frame((listPromoterCov)))
#list_promoter <-as.data.frame((listPromoterCov))
#write.csv(list_promoter,file=paste0("/usr/home/maloku/RiboProfiling_test/Results/","A1_list_promoter.csv"))
#png(filename="/usr/home/maloku/RiboProfiling_test/Results/trackPlotTSS_A1_prova.png",width=700, height=600)
#plot(list_promoter$X28_match.start,list_promoter$X28_match.values)
#lines(list_promoter$X28_match.start,list_promoter$X28_match.values)
#lines(list_promoter$X27_match.start,list_promoter$X27_match.values, col="red")
#lines(list_promoter$X29_match.start,list_promoter$X29_match.values, col="green")
#lines(list_promoter$X30_match.start,list_promoter$X30_match.values, col="blue")
#dev.off()
png(filename="/usr/home/maloku/RiboProfiling_test/bam_results/trackPlotTSS/trackPlotTSS.Rd_%03d_a1_prova.png",width=700, height=1000)
trackPlotTSS <-plotSummarizedCov(listPromoterCov)
trackPlotTSS
dev.off()

#keep only the match read sizes 27-33

alnGRanges <- alnGRanges[which(!is.na(match(alnGRanges$score,27:33)))]
countsData<-countShiftReads(exonGRanges[names(cdsPosTransc)], cdsPosTransc,
                                  alnGRanges,16,motifSize=3)
head(countsData[[1]],n=3)

#plot the coverage of the different features
png(filename="/usr/home/maloku/RiboProfiling_test/Results/countsPlot.Rd_%03d_A5_3.png", width=700, height=700)
listCountsPlots <- countsPlot(list(countsData[[1]]),grep("_counts$",colnames(countsData[[1]])),log2Bool=1)
listCountsPlots
dev.off()

### ribosome PAUSE scores on codon motifs
#for motifs of 3 nucleotides (1 codon)
listReadsCodon <- countsData[[2]]
dataset <-do.call(rbind, listReadsCodon)
write.csv(dataset,file=paste0("/usr/home/maloku/RiboProfiling_test/","A3_3_reads.csv"))

#for(i in length(listReadsCodon[1])){
  #print(listReadsCodon[i])
 # dat<-listReadsCodon[i]
  #print(dat)
  #write.xlsx(dat, file=paste0("/usr/home/maloku/RiboProfiling_test/bam_results/","a1_reads.xlsx"), sheetName = "dat", 
   #          col.names = TRUE, row.names = TRUE, append = FALSE)
  
#}

################## get the names of the expressed ORFs grouped by transcript

orfCoord <- cds[names(cds) %in% names(listReadsCodon)]

###### upload BSgenome.Athaliana.TAIR.TAIR9

genome <- BSgenome.Athaliana.TAIR.TAIR9
seqlevels(genome)
gSeq <- GenomeInfoDb::renameSeqlevels(genome,sub("Chr", "", GenomeInfoDb::seqlevels(genome)))
gSeq <- GenomeInfoDb::renameSeqlevels(gSeq,sub("M", "Mt", GenomeInfoDb::seqlevels(gSeq)))
gSeq <- GenomeInfoDb::renameSeqlevels(gSeq,sub("C", "Pt", GenomeInfoDb::seqlevels(gSeq)))
seqlevels(gSeq)

#codon frequency, coverage, and annotation

codonData <- codonInfo(listReadsCodon ,gSeq, orfCoord,motifSize=3)
data <- codonData
data <-as.data.frame(data)
write.csv(data, file=paste0("/usr/home/maloku/RiboProfiling_test/bam_results/","a9_codondata.csv"))
head(codonData[[1]], n=3)

##keep only genes with a minimum number of reads

codonUsage <- codonData[[1]]
codonCovMatrix <- codonData[[2]]
nbrReadsGene <- apply(codonCovMatrix, 1, sum)
ixExpGenes <- which(nbrReadsGene >= 50)
codonCovMatrix <- codonCovMatrix[ixExpGenes, ]

######Create CovMatrix , matrix of read density on codons #################

codonCovMatrixTransp <- t(codonCovMatrix)

#write.csv(codonCovMatrixTransp, file=paste0("/usr/home/maloku/RiboProfiling_test/bam_results/","a1cov_Matrix.csv"))
rownames(codonCovMatrixTransp) <- colnames(codonCovMatrix)
colnames(codonCovMatrixTransp) <- rownames(codonCovMatrix)
#write.csv(codonCovMatrixTransp, file=paste0("/usr/home/maloku/RiboProfiling_test/","A1_cov_Matrix.csv"))

###codons accumulating ribosome footprints-which are the most discriminant codon in the PCA analyses

png(filename="/usr/home/maloku/RiboProfiling_test/bam_results/pca/pca.Rd_%03d_a1_prova.png", width=700, height=700)
listPCACodonCoverage <- codonPCA(codonCovMatrixTransp, "codonCoverage")
print(listPCACodonCoverage[[2]])
#write.csv(listPCACodonCoverage[1], file=paste0("/usr/home/maloku/RiboProfiling_test/bam_results/","a23_listPCACodonCoverage.csv"))
dev.off()

#############PCA for amino acids which cause ribosome stalling##################

pca_all<-codonCovMatrixTransp
pca_all.pr<-prcomp(pca_all, center = TRUE, scale = TRUE, retx = TRUE)
summary(pcr_all.pr)
pca_all.pr$x

png(filename="/usr/home/maloku/RiboProfiling_test/Results//pca.Rd_%03d_A1_pca_prcomp_screenplot.png", width=700, height=700)
screeplot(pcr_all.pr, type = "l", npcs = 12, main = "plot of the first 12 PCs")
abline(h = 1, col="red", lty=5)
legend("topright", legend=c("Eigenvalue = 1"),
       col=c("red"), lty=5, cex=0.6)

cumpro <- cumsum(pcr_all.pr$sdev^2 / sum(pcr_all.pr$sdev^2))
plot(cumpro[0:12], xlab = "PC #", ylab = "Amount of explained variance", main = "Cumulative variance plot")
abline(v = 2, col="blue", lty=5)
abline(h = 0.88759, col="blue", lty=5)
legend("topleft", legend=c("Cut-off @ PC2"),
       col=c("blue"), lty=5, cex=0.6)
dev.off()

png(filename="/usr/home/maloku/RiboProfiling_test/Results/pca.k_means_A1_prova_screenplot.png", width=700, height=700)
wss <- (nrow(pca_all)-1)*sum(apply(pca_all,2,var))
for (i in 1:15) wss[i] <- sum(kmeans(pca_all,
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")
dev.off()

# From scree plot elbow occurs at k = 2
# Apply k-means with k=2

png(filename="/usr/home/maloku/RiboProfiling_test/Results/pca.Rd_%03d_A1_codonPCA.png", width=700, height=700)
pr_codon_all<-codonPCA(pca_all, "codonCoverage")
print(pr_codon_all)
dev.off()


#translate codons to AAs

aminoAcidsAAStringSet <- Biostrings::translate(DNAStringSet(colnames(codonUsage)))
aminoAcids <- unlist(strsplit(toString(aminoAcidsAAStringSet)[1], ", "))
#colnames(codonUsage[60])
ixaa <- grep("^L", aminoAcids)
#x<-(colnames(codonUsage))
#a<-grep("^TGT", x)
codonUsage_ixaa <- codonUsage[, ixaa]
codonCovMatrixTransp_ixaa <- codonCovMatrixTransp[ixaa, ]

#group by amino acids

aaX <- unique(aminoAcids[ixaa])

aminoAcidsUsage_ixaa <- matrix(nrow=nrow(codonUsage_ixaa), ncol=length(aaX))
aminoAcidsCov_ixaa <- matrix(ncol=ncol(codonCovMatrixTransp_ixaa), nrow=length(aaX))

for(ix in 1:length(aaX)){
  ixColsaa<-which(aminoAcids[ixaa] == aaX[ix])
  aminoAcidsUsage_ixaa[, ix] <- rowSums(codonUsage_ixaa[, ixColsaa])
  aminoAcidsCov_ixaa[ix, ] <- colSums(codonCovMatrixTransp_ixaa[ixColsaa, ])
}

colnames(aminoAcidsUsage_ixaa) <- aaX
rownames(aminoAcidsUsage_ixaa) <- rownames(codonUsage_ixaa)
colnames(aminoAcidsCov_ixaa) <- colnames(codonCovMatrixTransp_ixaa)
rownames(aminoAcidsCov_ixaa) <- aaX

codon_aa<-codonPCA(codonCovMatrixTransp_ixaa, "codonCoverage")
print(codon_aa)
aa<-prcomp(codon_aa$PCA_scores, center = TRUE, scale = TRUE, retx = TRUE)
print(aa)
