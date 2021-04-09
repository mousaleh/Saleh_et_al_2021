
#######################################################    ChEC-seq and ChIP-seq  plots ####
#                    ####
#                   #####
#########################


# Importing files ---------------------------------------------------------

ChEC<-read.csv("ChEC_ChIP-seq/All_ChEC_TSS_TES_log_changes.csv",
               header = T, stringsAsFactors = F)
ChIP<-read.csv("ChEC_ChIP-seq/All_X_ChIP_TSS_TES_log_changes.csv",
               header = T, stringsAsFactors = F)


# Statistical tests -------------------------------------------------------

# Make a dataframe to contain the p values
myp_val<-as.data.frame(grep("_16D", colnames(ChEC), value = T))
colnames(myp_val)[1]<-"Med_subunit"#changing the name of the column
myp_val$Med_subunit<-gsub("_16D", "",myp_val$Med_subunit) # removing _16D from the name
myp_val$Med16_up.pvalue<- 1 #create a column for ChEC signal at upregulated genes comparison (WT vs med16Del)
myp_val$Med16_down.pvalue<- 1 #create a column for ChEC signal at downregulated genes comparison (WT vs med16Del)

#loop to perform statistical test
for (i in seq(1,nrow(myp_val),1)){ #looping over the row number
  x<-paste(myp_val$Med_subunit[i],"16D", sep = "_") #fetching the subunit ChEC in med16 deletion
  y<-paste(myp_val$Med_subunit[i],"WT", sep = "_") #fetching the subunit ChEC in WT
  test.x<-wilcox.test(ChEC[c(which(ChEC$Cluster == "Med16_up")), grep(x,colnames(ChEC))], #conducting wilcoxon test and saing the result to a variable
                                          ChEC[c(which(ChEC$Cluster == "Med16_up")), grep(y,colnames(ChEC))])
  myp_val$Med16_up.pvalue[i]<-test.x$p.value # extracting the p-value only from the variable
  test.y<-wilcox.test(ChEC[c(which(ChEC$Cluster == "Med16_down")), grep(x,colnames(ChEC))],
                                          ChEC[c(which(ChEC$Cluster == "Med16_down")), grep(y,colnames(ChEC))])
  myp_val$Med16_down.pvalue[i]<-test.y$p.value
}


##pivoting the table
library(tidyverse)
ChEC_pivot<-gather(ChEC[,c(2:21)], key=Cluster, value="Log2_TSSR_TESR")
#changing the name of the column
colnames(ChEC_pivot)[2]<-"category"
ChEC_pivot$category<-factor(ChEC_pivot$category, levels = colnames(ChEC)[c(3:21)])


# Plotting Box-plots ------------------------------------------------------

#Plotting the box plots
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(ggpubr)

#plotting ChEc in WT and Med16 deletion
med_ChEC<-colnames(ChEC)[c(4:21)]
med_ChEC_p<-list()
temp_data<- ChEC_pivot %>% filter(category == "MED16_WT")
p1<-ggplot(temp_data, aes(Cluster, Log2_TSSR_TESR, fill = category))+geom_boxplot(size=0.25,outlier.size = 0.25)+theme_light()+scale_fill_brewer(palette="Paired")+scale_y_continuous(expand = c(0,0))+
  theme(panel.border = element_rect(colour = "black", size=1), legend.position = "none")+ylim(-5,15)+
  ylab("log2 (TSSR/TESR)")+xlab("Clusters")+ggtitle("Med16")
v<-1
for(i in seq(1,17,2)){
  temp_data<- ChEC_pivot %>% filter(category == med_ChEC[i] | category == med_ChEC[i+1])
  message(v)
  med_ChEC_p[[v]]<-local({
    v<-v
    ggplot(temp_data, aes(Cluster, Log2_TSSR_TESR, fill = category))+geom_boxplot(size=0.25,outlier.size = 0.25)+theme_light()+scale_fill_brewer(palette="Paired")+scale_y_continuous(expand = c(0,0))+
      theme(panel.border = element_rect(colour = "black", size=1), legend.position = "none")+ylim(-5,15)+
      ylab("log2 (TSSR/TESR)")+xlab("Clusters")+ggtitle(med_ChEC[i])
    #print(p1)
  })
  print(v)
  v<-v+1
  rm(temp_data)
}


#plotting in one figure
pdf("ChEC_ChIP-seq/All_ChEC_signals-test.pdf", width = 24, height = 18)
plot_grid(plotlist = med_ChEC_p,nrow = 2)
dev.off()

pdf("ChEC_ChIP-seq/Med16-ChEC_signals-test.pdf")
print(p1)
dev.off()
