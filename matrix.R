library(seqinr)

#Read file .fasta which contains genes and creates a list ''data''  

data <-read.fasta("Arabidopsis_thaliana.TAIR10.cdna.all.fa",as.string=TRUE, set.attributes=FALSE)

#It creates a vector with only gene names

genename <- attributes(data)

genenameonly <- unlist(genename, use.names =FALSE)

#It counts the lenght for each string in the list ''data'' and put  in the vector genelength

genelength <-rapply(data,nchar)

genelenghthonly <- as.vector(genelength)

#create a matrix

allmatrix <- cbind(genenameonly,genelenghthonly)

allmatrix_1 <- allmatrix[order(genenameonly),]

write.csv(allmatrix_1,file="/usr/home/maloku/Desktop/Davide_Ribo/tair10_genes_lengths.csv", quote=FALSE)
