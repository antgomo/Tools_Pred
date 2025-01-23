setwd("/Users/agomez/IMIDOMICS/cfDNA/")

#Boxplot
library(ggpubr)

annot<-read.delim("annot.csv",sep="\t")
dat2<-read.delim("cfDNA_results_batch2_new.csv",sep="\t",head=F)
#dat2<-read.csv("Results_5cpgs_batch2.csv",head=F)

#dat2<-dat2[!(dat2$V3==""),]

colnames(dat2)[1]<-"multiplex_index"
dat2<-merge(dat2,annot, by="multiplex_index")
colnames(dat2)[2]<-"cfDNA"

#colnames(dat2)<-c("V1","cfDNA","disease")

p <- ggboxplot(dat2[-4,c("cfDNA","disease")], x = "disease", y = "cfDNA",##change depending on set of genes
               color = "disease",
               add = "jitter",
               legend="")

y.df<-dat2[-4,c("cfDNA","disease")]
result<-aov(y.df[,1]~as.factor(y.df$disease))

res<-summary(result)
pval<-res[[1]][["Pr(>F)"]][1]
pval=0.03537505
pval<-signif(pval, digits=3)

p<-p + xlab("") + ggtitle(paste("Batch2 cfDNA quantification",pval,sep=" "))+
  rotate_x_text(angle = 45)+
  geom_hline(yintercept = median(dat2[dat2$disease=="RA","cfDNA"]), linetype = 2,color="lightblue") # Add horizontal line at bas
# Specify the comparisons you want
my_comparisons <- list(  c("CTRL", "RA"))

p<- p + stat_compare_means(method="t.test",comparisons = my_comparisons,
                           #
                           label = "p.signif",###only signifficance symbol stars
                           # method.args = list(var.equal = TRUE),##make t-test similar to ANOVA 2 groups,
                           symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                               symbols = c("****", "***", "**", "*", "ns"))
)
p



###corrrelation batch1

#dat1<-read.delim("cfDNA_results_B1_GEM.csv",head=F,sep="\t")
dat1<-read.csv("cfDNA_results_B1_GEM.csv",head=F)

##select dilutions
dat1<-dat1[c(2,3,7,9,11),]

##order

my.data<-dat1[c(1,5,4,3,2),1:2]
my.data$V3<-c(0.1,0.5,1,10,100)  
#my.data<-my.data[-5,]
  
  mod1 <-lm(as.numeric(my.data[,1])~as.numeric(my.data[,3]))## 
  modsum<-summary(mod1)
  
  r2 <- cor.test(as.numeric(my.data[,1]),as.numeric(my.data[,3]),method="spearman")$estimate
  my.p<-cor.test(as.numeric(my.data[,1]),as.numeric(my.data[,3]),method="spearman")$p.value
  my.p<-signif(my.p, digits=3)
  
  
  mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
  plot(as.numeric(my.data[,1])~as.numeric(my.data[,3]), main=paste("Correlation Dilutions GEM alignment",sep=" "), 
       xlab="Dilutions", ylab="cfDNA Unmeth Percentage",
       pch=20)
  legend('topleft', legend = mylabel, bty = 'n')
  abline(mod1, col="red")
  

