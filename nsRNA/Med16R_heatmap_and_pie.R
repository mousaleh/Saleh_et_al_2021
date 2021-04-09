#importing text files
med16.AID<-read.delim("/Users/zentlab/Desktop/Thesis project/paper_1_Github/Data/DEGs/med16AID.txt", header = T, stringsAsFactors = F)
med15.AID<-read.delim("/Users/zentlab/Desktop/Thesis project/paper_1_Github/Data/DEGs/Med15_AID.txt", header = T, stringsAsFactors = F)
med15_16.AID<-read.delim("/Users/zentlab/Desktop/Thesis project/paper_1_Github/Data/DEGs/Med15_16_AID.txt", header = T, stringsAsFactors = F)
med16_mot1.AID<-read.delim("/Users/zentlab/Desktop/Thesis project/paper_1_Github/Data/DEGs/Med16_Mot1_AID.txt", header = T, stringsAsFactors = F)
med16.Del<-read.delim("/Users/zentlab/Desktop/Thesis project/paper_1_Github/Data/DEGs/med16DvsWt.txt", header = T, stringsAsFactors = F)

# making the Med16 clusters list ------------------------------------------

#Making a list of genes significantly altered in Med16-AID and/or Med16-Del
Med16_genes<-as.character(med16.AID$Transcript.ID[c(which(med16.AID$padj < 0.05))])
Med16_genes<-append(Med16_genes,as.character(med16.Del$Transcript.ID[c(which(med16.Del$padj < 0.05))]))

#remove duplicated genes
Med16_genes<-Med16_genes[!duplicated(Med16_genes)]

#remove three genes that span or related to Med16 deletion Leu2 (YCL018W), Med16 (YNL236W), YNL234W
Med16_genes<-Med16_genes[-c(grep("YCL018W|YNL236W|YNL234W", Med16_genes))]

#adding a prefix to colnames
colnames(med16.AID)[c(1:6)]<-paste("med16.AID",colnames(med16.AID)[c(1:6)], sep = "_")
colnames(med16.Del)[c(1:6)]<-paste("med16.Del",colnames(med16.Del)[c(1:6)], sep = "_")

#one dataframe for both conditions
All.med16<-merge(med16.AID[,c(2,7)], med16.Del[,c(2,7)], by="Transcript.ID")
All.med16<-merge(as.data.frame(Med16_genes), All.med16, by.x="Med16_genes", by.y="Transcript.ID", all.x=TRUE)

#retain Genes with concordant effects in AID and del vs WT
All.med16 <- All.med16 [(All.med16$med16.AID_log2FoldChange > 0 & All.med16$med16.Del_log2FoldChange > 0 ) | (All.med16$med16.AID_log2FoldChange < 0 & All.med16$med16.Del_log2FoldChange < 0 ),]

#making a column with the gene categories and adding it to the dataframe
All.med16$Category<-"Med16_up"
All.med16$Category[c(which(rowMeans(All.med16[,c(2,3)]) < 0))]<-"Med16_down"


# Heatmap of Med16/AID-Del at Med16-gene clusters -------------------------

library(pheatmap)
library(RColorBrewer)
m2<-as.matrix(All.med16[c(order(All.med16$Category)),c(2,3)])
palette_length <- 1000
my_pal <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(palette_length))
my_breaks <- c(seq(min(m2), 0, length.out = ceiling(palette_length / 2) + 1), 
               seq(max(m2) / palette_length, 1, length.out = floor(palette_length / 2)))
pdf("/Users/zentlab/Desktop/Thesis project/paper_1_Github/Med16-filtered-heatmap.pdf")
pheatmap(m2, cluster_rows = F, cluster_cols = F, breaks = my_breaks, color = my_pal, show_rownames = F)
dev.off()


# Heatmap of vORF ---------------------------------------------------------

#changing the name of log fold change column
colnames(med15.AID)[2]<-"Med15.AID"
colnames(med15_16.AID)[2]<-"Med15_16.AID"
colnames(med16_mot1.AID)[2]<-"Med16_Mot1.AID"

##making one dataframe with all DEGs
myDEGs<-merge(med16.AID[,c(2,7)],merge(med15.AID[,c(2,7)], merge(med15_16.AID[,c(2,7)],med16_mot1.AID[,c(2,7)])))
colnames(myDEGs)<-gsub("_log2FoldChange", "", colnames(myDEGs))# making long column names shorter

#importing the vORF
vORF<-read.csv("verified_ORFs.csv", header = F, stringsAsFactors = F)
vORF.DEGs<-merge(vORF, myDEGs, all.x = T, by.x = "V2", by.y = "Transcript.ID")

#plotting heatmaps

Med15.h<-as.matrix(vORF.DEGs[,c(6,7,8)])
Med15.h<-na.omit(Med15.h)
Mot1.h<-as.matrix(vORF.DEGs[,c(6,9)])
Mot1.h<-na.omit(Mot1.h)


library(pheatmap)
library(RColorBrewer)

#Med15
m2<-Med15.h
palette_length <- 1000
my_pal <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(palette_length))
my_breaks <- c(seq(min(m2), 0, length.out = ceiling(palette_length / 2) + 1), 
               seq(max(m2) / palette_length, 1, length.out = floor(palette_length / 2)))
pdf("Med15-vORF-heatmap_clustered.pdf")
pheatmap(m2, cluster_rows = T, cluster_cols = F, breaks = my_breaks, color = my_pal, show_rownames = F)
dev.off()


##Mot1
m2<-Mot1.h
palette_length <- 1000
my_pal <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(palette_length))
my_breaks <- c(seq(min(m2), 0, length.out = ceiling(palette_length / 2) + 1), 
               seq(max(m2) / palette_length, 1, length.out = floor(palette_length / 2)))
pdf("Mot1-vORF-heatmap_clustered.pdf")
pheatmap(m2, cluster_rows = T, cluster_cols = F, breaks = my_breaks, color = my_pal, show_rownames = F)
dev.off()



# Pie-chart with Hahn categories ------------------------------------------


#importing Steve Hahn Categories
library('tidyverse')
library(viridis)
hahn20<-read.delim("/Users/zentlab/Desktop/Thesis project/paper_1_Github/Donczew-et-al-elife-2020-categories.txt",
                   header = T, stringsAsFactors = F)
#new Hahn categories enrichment
colnames(hahn20)[1]<-"Transcript.ID"
#annotating  genes by new Donczew et al dataset
Category.genes1<-merge(All.med16[,c(1,4)], hahn20, by.x = "Med16_genes", by.y = "Transcript.ID", all.x=TRUE)
#making a bar chart of the new categories
mycounts<-Category.genes1[,c(2,4)]
mycounts[,2]<-as.character(mycounts[,2])
#removing empty and NA rows and replacing them by "no call"
mycounts[(is.na(mycounts$class)| mycounts$class == ""),2]<-"no call"
#converting the different categories to factors
mycounts[,2]<-factor((mycounts[,2]), levels = c("no call", "CR", "TFIID"))
#splitting the data into Med16_up and MEd16_down data
data.up <- mycounts[c(which(mycounts$Category == "Med16_up")),]
data.down <- mycounts[c(which(mycounts$Category == "Med16_down")),]
#summarizing each by the different classes
data.up <- data.up %>% group_by(Category,class) %>% summarize(n=n())
data.down <- data.down %>% group_by(Category,class) %>% summarize(n=n())

#plotting the pie-chart for Med16_up
Up_pie<-ggplot(data.up, aes(x="", y=n,fill = class, group = class)) +
  geom_bar(stat="identity", width=1)+ coord_polar("y", start=0) + theme_void()+scale_fill_manual(values = viridis(4))+
  geom_text(aes(x="", y = n, label=n))
pdf("Hahn20-pie-up1.pdf")
print(Up_pie)
dev.off()

#plotting the pie-chart for Med16_down
down_pie<-  ggplot(data.down, aes(x="",y=n, fill = class, group = class)) +
  geom_bar(stat="identity", width=1)+ coord_polar("y", start=0) + theme_void()+scale_fill_manual(values = viridis(4))+
  geom_text(aes(x="", y = n, label=n))
pdf("Hahn20-pie-down1.pdf")
print(down_pie)
dev.off()

# Med16 gene clusters box-plots -------------------------------------------

##clusters boxplots

library(tidyr)
myDEGs<-merge(Med16_genes, myDEGs, all.x=T, by.x="Med16_genes", by.y="Transcript.ID")
#pivoting with clusters
all.AID.pivot<-gather(myDEGs[,c(2,4:7)], key = "Category")
colnames(all.AID.pivot)[1]<-"Cluster"
all.AID.pivot$Category<-factor(all.AID.pivot$Category, levels = unique(all.AID.pivot$Category)) #ensuring the order of categories in plots

#plotting
library(ggplot2)
pdf("Med15-boxplots.pdf")
ggplot(all.AID.pivot[c(which(all.AID.pivot$Category == "med16.AID"|all.AID.pivot$Category == "Med15.AID"|all.AID.pivot$Category == "Med15_16.AID")),], aes(Cluster, value, fill = Category))+geom_boxplot(size=0.25,outlier.size = 0.25)+theme_light()+scale_y_continuous(expand = c(0,0))+
  theme(panel.border = element_rect(colour = "black", size=1), axis.text.x = element_text(face = "bold", color = "black"))+
  ylab("Log2 Fold Change")
dev.off()


pdf("Mot1-boxplots.pdf")
ggplot(all.AID.pivot[c(which(all.AID.pivot$Category == "med16.AID"|all.AID.pivot$Category == "Med16_Mot1.AID")),], aes(Cluster, value, fill = Category))+geom_boxplot(size=0.25,outlier.size = 0.25)+theme_light()+scale_y_continuous(expand = c(0,0))+
  theme(panel.border = element_rect(colour = "black", size=1), axis.text.x = element_text(face = "bold", color = "black"))+
  ylab("Log2 Fold Change")
dev.off()


# Box-plots Statistical test ----------------------------------------------

#making a dataframe with only Med15, Med16 and Med15/Med16-AID values
med15.box <- all.AID.pivot[c(which(all.AID.pivot$Category == "med16.AID"|all.AID.pivot$Category == "Med15.AID"|all.AID.pivot$Category == "Med15_16.AID")),]

#performing statistical tests with multiple tests correction

attach(med15.box[c(which(med15.box$Cluster == 1)),])
pairwise.wilcox.test(value, Category, p.adjust.method = "holm",
                     paired = FALSE)
detach()

attach(med15.box[c(which(med15.box$Cluster == 2)),])
pairwise.wilcox.test(value, Category, p.adjust.method = "holm",
                     paired = FALSE)
detach()

#statistical tests on Med16/Mot1-AID
wilcox.test(myDEGs$med16.AID[c(which(myDEGs$Category == "Med16_down" ))],myDEGs$Med16_Mot1.AID[c(which(myDEGs$Category == "Med16_down" ))])
wilcox.test(myDEGs$med16.AID[c(which(myDEGs$Category == "Med16_up" ))],myDEGs$Med16_Mot1.AID[c(which(myDEGs$Category == "Med16_up" ))])


# vORFs heatmap -----------------------------------------------------------


