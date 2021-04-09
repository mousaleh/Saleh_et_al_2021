
#######################################################    nsRNA analysis  ####
                                                      #                    ####
                                                      #                   #####
                                                      #########################


library(DESeq2)
library(ggplot2)
library(dplyr)

##getting raw counts
Budding_count<-read.csv("/Users/zentlab/Desktop/Thesis project/paper_1_Github/Budding_raw_count1.csv", header = T)
Pombe_count<-read.csv("/Users/zentlab/Desktop/Thesis project/paper_1_Github/Pombe_raw_count1.csv", header = T)


# Preparing DEseq2 object -------------------------------------------------


#preparing count matrix
cts <- Budding_count
rownames(cts) <- cts[,1]
cts[,1] <- NULL
cts <- round(cts,0)
cts <- as.matrix(cts)
coldata<-as.data.frame(colnames(cts))
coldata$type<-'paired-end'
coldata$genotype_treatment<-"GT"
coldata$genotype_treatment[c(1:12)] <- rep(rep(c("Med16.AID_IAA","Med16.AID_DMSO","Wt","med16.Del"),each=3),1)
coldata$genotype_treatment[c(13:18)] <- rep(rep(c("Med15_16.AID_DMSO","Med15_16.AID_IAA"),each=1),3)
coldata$genotype_treatment[c(19:24)] <- rep(rep(c("Med16_Mot1.AID_DMSO","Med16_Mot1.AID_IAA"),each=1),3)
coldata$genotype_treatment[c(25:30)] <- rep(rep(c("Med15.AID_DMSO","Med15.AID_IAA"),each=1),3)
coldata$genotype_treatment<-factor(coldata$genotype_treatment, levels = unique(coldata$genotype))
colnames(coldata)[1]<-"col.data"


# Spike-in Normalization --------------------------------------------------


#getting spike-in data
cts1 <- Pombe_count
rownames(cts1) <- cts1[,1]
cts1[,1] <- NULL
cts1 <- round(cts1,0)
cts1 <- as.matrix(cts1)

#getting size factors from spike-in count Matrix
x<-estimateSizeFactorsForMatrix(counts=cts1)

#creating the dESeq2 object
dds<- DESeqDataSetFromMatrix(countData = cts,
                             colData = coldata,
                             design = ~genotype_treatment)

#using Spike-in size factors to normalize library
sizeFactors(dds) <- x


# PCA ---------------------------------------------------------------------



##calculating PCA components and plotting PCA

#applying the variance stabilization function 
vsd <- vst(dds, blind=FALSE)


#Med16_PCA

pdf('/Users/zentlab/Desktop/Thesis project/paper_1_Github/Med16_PCA2.pdf')
pcaData<-plotPCA(vsd[,c(1:12)], intgroup="genotype_treatment", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=genotype_treatment)) +
  geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+ theme_minimal()+ 
  coord_fixed()
dev.off()

#Med15_PCA

pdf('/Users/zentlab/Desktop/Thesis project/paper_1_Github/Med15_PCA2.pdf')
pcaData<-plotPCA(vsd[,c(13:18,25:30)], intgroup="genotype_treatment", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=genotype_treatment)) +
  geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + theme_minimal()+
  coord_fixed()
dev.off()

#Med16/Mot1-AID PCA

pdf('/Users/zentlab/Desktop/Thesis project/paper_1_Github/Med16_Mot1-AID_PCA2.pdf')
pcaData<-plotPCA(vsd[,c(19:24)], intgroup="genotype_treatment", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=genotype_treatment)) +
  geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + theme_minimal()+
  coord_fixed()
dev.off()

