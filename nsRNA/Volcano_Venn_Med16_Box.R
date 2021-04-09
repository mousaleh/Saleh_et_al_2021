#Importing Differential expression genes
med16.AID<-read.delim("/Users/zentlab/Desktop/Thesis project/paper_1_Github/Data/DEGs/med16AID.txt", header = T, stringsAsFactors = F)
med15.AID<-read.delim("/Users/zentlab/Desktop/Thesis project/paper_1_Github/Data/DEGs/Med15_AID.txt", header = T, stringsAsFactors = F)
med15_16.AID<-read.delim("/Users/zentlab/Desktop/Thesis project/paper_1_Github/Data/DEGs/Med15_16_AID.txt", header = T, stringsAsFactors = F)
med16_mot1.AID<-read.delim("/Users/zentlab/Desktop/Thesis project/paper_1_Github/Data/DEGs/Med16_Mot1_AID.txt", header = T, stringsAsFactors = F)
med16.Del<-read.delim("/Users/zentlab/Desktop/Thesis project/paper_1_Github/Data/DEGs/med16DvsWt.txt", header = T, stringsAsFactors = F)


# Volcano plots -----------------------------------------------------------


#making the volcano plots

#med16_AID volcano plots
med16.AID<-med16.AID[c(which(med16.AID$padj != 0)),]
med16.AID$log10.p <- -log10(med16.AID$padj)
med16.AID$significance<-"NS"
med16.AID$significance[(med16.AID$padj < 0.05)] <- "FDR"
med16.AID$significance <- factor(med16.AID$significance, levels=c("NS", "FDR"))
volc_med16.AID = ggplot(med16.AID, aes(log2FoldChange, log10.p, group_by(significance))) + #volcanoplot with log2Foldchange versus pvalue
  geom_point(aes(col=significance), alpha=1, size = 1) + #add points colored by significance
  #choose specific shape to certain points manually 
  scale_shape_manual(values=c(1,1)) + labs(title="Med16-AID (IAA vs DMSO)")+
  #Choose which colours to use; otherwise, ggplot2 choose automatically (order depends on how factors are ordered in toptable$Significance)
  scale_color_manual(values=c(NS="grey30", FDR="red1"), labels=c(NS="NS", FDR=paste("FDR adj-p <", 0.05, sep=""))) +
  theme_classic() + xlim(-5.4, 5.1)+
  #Add a horizontal line for P-value cut-off
  geom_hline(yintercept=(-log10(0.05)), linetype = 'longdash', color = 'black', size=0.4) 
pdf("/Users/zentlab/Desktop/Thesis project/paper_1_Github/med16.AID1.pdf")
print(volc_med16.AID)
dev.off()

#med16_deletion vs Wt volcano plots
med16.Del<-med16.Del[c(which(med16.Del$padj != 0)),]
med16.Del$log10.p <- -log10(med16.Del$padj)
med16.Del$significance<-"NS"
med16.Del$significance[(med16.Del$padj < 0.05)] <- "FDR"
med16.Del$significance <- factor(med16.Del$significance, levels=c("NS", "FDR"))
volc_med16.Del = ggplot(med16.Del, aes(log2FoldChange, log10.p, group_by(significance))) + #volcanoplot with log2Foldchange versus pvalue
  geom_point(aes(col=significance), alpha=1, size = 1) + #add points colored by significance
  #choose specific shape to certain points manually 
  scale_shape_manual(values=c(1,1)) + labs(title="Med16-AID (IAA vs DMSO)")+
  #Choose which colours to use; otherwise, ggplot2 choose automatically (order depends on how factors are ordered in toptable$Significance)
  scale_color_manual(values=c(NS="grey30", FDR="red1"), labels=c(NS="NS", FDR=paste("FDR adj-p <", 0.05, sep=""))) +
  theme_classic() + xlim(-5.4, 5.1)+
  #Add a horizontal line for P-value cut-off
  geom_hline(yintercept=(-log10(0.05)), linetype = 'longdash', color = 'black', size=0.4) 
pdf("/Users/zentlab/Desktop/Thesis project/paper_1_Github/med16.Del1.pdf")
print(volc_med16.Del)
dev.off()


# Venn Diagrams -----------------------------------------------------------

#make a list of gene identifiers of each category
med16.up<-med16.AID$Transcript.ID[c(which(med16.AID$log2FoldChange > 0 & med16.AID$significance == "FDR"))]
med16D.up<-med16.Del$Transcript.ID[c(which(med16.Del$log2FoldChange > 0 & med16.Del$significance == "FDR"))]
med15.up<-med15.AID$Transcript.ID[c(which(med15.AID$log2FoldChange > 0))]
med15_16.up<-med15_16.AID$Transcript.ID[c(which(med15_16.AID$log2FoldChange > 0))]

med16.dwn<-med16.AID$Transcript.ID[c(which(med16.AID$log2FoldChange < 0 & med16.AID$significance == "FDR"))]
med15.dwn<-med15.AID$Transcript.ID[c(which(med15.AID$log2FoldChange < 0))]
med15_16.dwn<-med15_16.AID$Transcript.ID[c(which(med15_16.AID$log2FoldChange < 0))]
med16D.dwn<-med16.Del$Transcript.ID[c(which(med16.Del$log2FoldChange < 0 & med16.Del$significance == "FDR"))]

library(VennDiagram)
library(viridis)
#initial venns for feeding into venneuler

#med16 (upregulated)
venn.diagram(
  x = list(med16D.up, med16.up),
  category.names = c("med16D" , "Med16_AID"),
  filename = '/Users/zentlab/Desktop/Thesis project/mediator images/med16Del_med16-AID_up.tiff',
  output=TRUE,resolution =500,
  imagetype = "tiff", units = "px", compression =
    "lzw",
  output=TRUE,
  col=viridis(2),
  fill = viridis(2)
)

#med16 (downregulated)
venn.diagram(
  x = list(med16D.dwn, med16.dwn),
  category.names = c("med16D" , "Med16_AID"),
  filename = '/Users/zentlab/Desktop/Thesis project/mediator images/med16Del_med16-AID_dwn.tiff',
  output=TRUE,resolution =500,
  imagetype = "tiff", units = "px", compression =
    "lzw",
  output=TRUE,
  col=viridis(2),
  fill = viridis(2)
)


##med15 and Med16 concordance
venn.diagram(
  x = list(med16.dwn, med15.dwn),
  category.names = c("Med16-AID" , "Med15_AID"),
  filename = '/Users/zentlab/Desktop/Thesis project/mediator images/med15_16-AID_dwn.tiff',
  output=TRUE,resolution =500,
  imagetype = "tiff", units = "px", compression =
    "lzw",
  output=TRUE,
  col=viridis(2),
  fill = viridis(2)
)

#final figure Venns

library(venneuler)
MyVenn <- venneuler(c(Med16D_up=578,Med16AID_up=340, 
                      "Med16D_up&Med16AID_up"=100))
MyVenn$labels <- c("med16D\n478","Med16-AID\n240")
plot(MyVenn)
text(0.475,0.51,"100")

MyVenn <- venneuler(c(Med16D_dwn=246,Med16AID_dwn=44, 
                      "Med16D_dwn&Med16AID_dwn"=34))
MyVenn$labels <- c("med16D\n212","Med16-AID\n10")
plot(MyVenn)
text(0.475,0.51,"34")

MyVenn <- venneuler(c(Med16AID_up=340,Med15_up=297,
                      "Med16AID_up&Med15_up"=32))
MyVenn$labels <- c("Med16-AID","Med15-AID")
plot(MyVenn)
text(0.475,0.51)

MyVenn <- venneuler(c(Med16AID_dwn=44,Med15_dwn=183,
                      "Med16AID_dwn&Med15_dwn"=10))
MyVenn$labels <- c("Med16-AID","Med15-AID")
plot(MyVenn)
text(0.475,0.51)

# med16 DEGS box-plots ----------------------------------------------------
library(ggplot2)

#making a dataframe with both Med16-AID and Med16.Del
med16.all <- merge(med16.AID[,c(2,7)], med16.Del[,c(2,7)], by = "Transcript.ID")

#renaming the columns
colnames(med16.all)[c(2:3)]<-c("med16.AID", "med16.Del")
#use the list of genes DE in Med16_AID or Med16-deletion to make a dataframe with category labelled
med16.DE<-list("med16.AID_up"= med16.up, "med16.AID_down"= med16.dwn, "med16.DEL_up"= med16D.up, "med16.DEL_down"= med16D.dwn)

#making a loop to create a dataframe that has the category and names of genes
for (l in seq(length(med16.DE))){
  if (!exists("mydf")){#checking if mydf elists to add to it or create it
    mydf<-as.data.frame(med16.DE[[l]]) #making a list into a dataframe
    mydf$category<-names(med16.DE[l]) # use the list name as a string in the category column
    colnames(mydf)[1]<-"Transcript.ID"
  }
  if (exists("mydf")){
    y<-as.data.frame(med16.DE[[l]])
    y$category<-names(med16.DE[l])
    colnames(y)[1]<-"Transcript.ID"
    mydf<-rbind(mydf,y) #combining the new list to the mydf dataframe
  }
  }

#adding differential expression to the dataframe
Med16.DEGs<-merge(mydf, med16.all, all.x = T)

##pivoting the dataframe
library(tidyverse)
Med16.all.pivoted<-gather(Med16.DEGs[,c(2:4)], key = "category")
colnames(Med16.all.pivoted)[2]<-"Experiment"

#plotting
#Med16_DEl DE

pdf("/Users/zentlab/Desktop/Thesis project/paper_1_Github/Med16_Del_sig_genes-Med16_nsRNA-boxplot.pdf")
ggplot(Med16.all.pivoted[c(which(Med16.all.pivoted$category == "med16.DEL_up" |Med16.all.pivoted$category == "med16.DEL_down")),],
       aes(category, value, fill = Experiment))+geom_boxplot(size=0.25,outlier.size = 0.25)+theme_light()+
  theme(panel.border = element_rect(colour = "black", size=1), axis.text.x = element_text(face = "bold", color = "black"))+
  ylab("Log2 Fold Change")
dev.off()

#Med16_AID DE

pdf("/Users/zentlab/Desktop/Thesis project/paper_1_Github/Med16_AID_sig_genes-Med16_nsRNA-boxplot.pdf")
ggplot(Med16.all.pivoted[c(which(Med16.all.pivoted$category == "med16.AID_up" |Med16.all.pivoted$category == "med16.AID_down")),],
       aes(category, value, fill = Experiment))+geom_boxplot(size=0.25,outlier.size = 0.25)+theme_light()+
  theme(panel.border = element_rect(colour = "black", size=1), axis.text.x = element_text(face = "bold", color = "black"))+
  ylab("Log2 Fold Change")
dev.off()