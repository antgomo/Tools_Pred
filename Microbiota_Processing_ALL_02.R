library(data.table)
library(tidyselect)

setwd("/Users/agomez/IMIDOMICS/Leakage/")


##ctrls
load("Microbiota_CTRL_ALL.RData")

###IMID
load("Microbiota_P6_Gut_ALL.RData")


###merge both

all_data<-c(lst,ctrls_data)

######now, we have to check if both reads are mapping against the same organism
####AVOID this step!!!

###01 Reads to 01
###02 Filter to more than 85 length
###03 Consider singletons
###not clear, better results in 1 million filtering whithout singletons


filt<-all_data
rm(all_data)
gc()
results.list<-list()

for (i in 1:length(filt)){ 
  
  ids.list<-strsplit(filt[[i]]$qseqid,"/")
  
  ids<-sapply(ids.list, function(x) return(x[1]))
  
  filt[[i]]$ids<-ids
  
  ##split df
  
  read1<-filt[[i]][grep("/1",filt[[i]]$qseqid),c("qseqid","sseqid","ids")]
  read2<-filt[[i]][grep("/2",filt[[i]]$qseqid),c("qseqid","sseqid","ids")]
  
  ##merge each read by common name to get both mates in same line
  big.df<-merge(read1[,c("qseqid","sseqid","ids")],read2[,c("qseqid","sseqid","ids")],by="ids",allow.cartesian=TRUE)
  
  ######here we can check for 1
  ###check if both mates have mapped to same organism
  big.df$CHECK<-ifelse(big.df$sseqid.x==big.df$sseqid.y,"TRUE","FALSE")
  ####get only unambiogus matches aka. both mates same organism
  true.matches<-big.df[big.df$CHECK=="TRUE",]
  ##collect for each sample
  results.list[[i]]<-as.data.frame(table(true.matches$sseqid.x))
  
}


###########try to plot abundance

library(dplyr)
library(purrr)
####get all results in one matrix
#res.df<-results.list[-30] %>% purrr::reduce(left_join, by = "Var1")##for PS
res.df<-results.list %>% purrr::reduce(left_join, by = "Var1")##for UC

#colnames(res.df)[-1]<-names(filt)[-30]##for PS
colnames(res.df)[-1]<-names(filt)##for UC

colnames(res.df)[1]<-"V1"
rownames(res.df)<-res.df$V1
res.df<-res.df[,-1]


###load species header

species<-read.delim("Gut_species1.txt",sep=" ", header=F)
species<-species[,1:3]
species$ID<-paste(species$V2,species$V3,sep="_")


##################

###order using samplesheet
####load samplesheet
samplesheet<-read.delim("../classes/SampleTable_FINAL.csv")
samplesheet<-samplesheet[samplesheet$GSL.ID %in% colnames(res.df),]
res.df<-res.df[,match(samplesheet$GSL.ID,colnames(res.df))]

######change NAs to 0
res.df[is.na(res.df)] <- 0

save(res.df,file="Matrix_counts_Microbiota_RA_CTRLs.RData")


#load counts

###try to normalize mapped reads first

######consider leakage as the total of reads mapped in gut microbiome
leak<-as.data.frame(colSums(res.df))

leak$GSL.ID<-rownames(leak)
colnames(leak)[1]<-"Leak"

#counts<-merge(all.counts,leak, by="GSL.ID")
###order using samplesheet
samplesheet<-read.delim("../classes/SampleTable_FINAL.csv")

samplesheet<-samplesheet[samplesheet$GSL.ID %in% leak$GSL.ID,]
leak<-leak[match(samplesheet$GSL.ID,as.character(leak$GSL.ID)),]

####use 1000000
leak$LEAK_PERC<-(leak$Leak/1000000)*100

my.df<-merge(samplesheet[,c("GSL.ID","imid","activity")],leak,by="GSL.ID")
my.df$ACT<-ifelse(is.na(my.df$activity),"CTRL",ifelse(my.df$activity=="HI","HI","LO"))

my.df$ACT <- factor(my.df$ACT,levels = c("HI", "LO", "CTRL"))

save(my.df,file="Leakage_RA_CTRLS.RData")

##try to normalize
library(BBmisc)
#“standardize”: Center and scale
my.df$LEAK_PERC<-BBmisc::normalize(my.df$LEAK_PERC, method = "standardize", range = c(0, 1))
###Boxplot
library(ggpubr)

####remove ctrl outliers

my.df<-my.df[my.df$LEAK_PERC<2,]###UC samples

p <- ggboxplot(my.df, x = "ACT", y = "LEAK_PERC",##change depending on set of genes
               color = "imid",
               add = "jitter",
               legend="",
               outlier.shape = ""
               #ylim=c(-0.4,2)
)



p<-p+ ggtitle("Leakage Activity RA CTRLs")+
  rotate_x_text(angle = 45)
# Specify the comparisons you want
my_comparisons <- list( c("HI", "LO"),c("HI","CTRL"),c("LO","CTRL"))
##add covariates

samplesheet$ACT<-ifelse(is.na(samplesheet$activity),"CTRL",ifelse(samplesheet$activity=="HI","HI","LO"))
samplesheet$batch<-paste("B",samplesheet$new.plate,sep="")
samplesheet$treat<-ifelse(samplesheet$activity=="HI","HI",ifelse(samplesheet$activity=="LO","LO","CTRL"))
samplesheet$treat[is.na(samplesheet$treat)]<-"CTRL"
samplesheet$Sex<-as.factor(ifelse(samplesheet$sex=="FEMALE","F","M"))
treat<-as.character(samplesheet$treat)
pheno<-samplesheet[,c("ACT","GSL.ID","Sex","age","batch")]
rownames(pheno)<-pheno$GSL.ID

my.df<-merge(pheno[,-1],my.df,by="GSL.ID")
###normalize alpha Inverse normal transformation
#Alpha.df$y <- qnorm((rank(Alpha.df$inverse_simpson, na.last="keep") - 0.5) / sum(!is.na(Alpha.df$inverse_simpson)))
#y<-as.data.frame(y)




#####to put on stats
summary(lm(LEAK_PERC ~ relevel(as.factor(my.df$ACT),"CTRL") +Sex + age + batch , data = my.df))





p<- p + stat_compare_means(method="wilcox.test",comparisons = my_comparisons,
                           #
                
                           #label = "p.signif",###only signifficance symbol stars
                           method.args = list(var.equal = T),##make t-test simila to ANOVA 2 groups,
                           
)
p
#####add linear model pvals instead wilcox
stat.test <- compare_means(
  LEAK_PERC ~ ACT, data = my.df,
  method = "wilcox.test"
)

stat.test$p<-c(0.89,0.025232 ,0.108679)

stat.test <- stat.test %>%
  mutate(y.position = c(2, 4, 3))

p + stat_pvalue_manual(stat.test, label = "p")




write.csv(my.df,"Leakage_Score_UC_CTRL_Fixed.csv")
###################
####GLOBAL MEASURES
#####


#Microbial diversity within a sample was determined using the richness and alpha diversity
#indices. Richness was defined as the total number of distinct taxa in a sample. We use Inverse
#Simpson's formula incorporating richness and evenness components to compute alpha diversity


#To measure sample to-sample dissimilarities between microbial communities we use Bray-Curtis beta diversity
#index accounting for both changes in the abundances of the shared taxa and account

#To test for differences in alpha diversity between disease groups, we fit the following analysis of
#covariance (ancova) model
#𝑎𝑙𝑝ℎ𝑎_𝑛𝑜𝑟𝑚 ~ 𝑆𝑒𝑥 + 𝐴𝑔𝑒 + 𝑇𝑒𝑐ℎ𝑛𝑖𝑐𝑎𝑙 𝑐𝑜𝑣𝑎𝑟𝑖𝑎𝑡𝑒𝑠 + 𝐷𝑖𝑠𝑒𝑎𝑠𝑒 𝑠𝑡𝑎𝑡𝑢𝑠

##Technical covariates include: RIN, Batch (Plate_number), Concentration



##########Diversity
library(vegan)

##Simpson: The probability that two randomly chosen individuals are the same species.
##Inverse Simpson: This is a bit confusing to think about. 
  #Assuming a theoretically community where all species were equally abundant, this would be the number of species needed to have the same Simpson index value 
  #for the community being analyzed.


library(microbiome)

tab <- microbiome::diversity(res.df, index = "all")


Alpha.df<-merge(pheno,tab,by="row.names")
###normalize alpha Inverse normal transformation
#Alpha.df$y <- qnorm((rank(Alpha.df$inverse_simpson, na.last="keep") - 0.5) / sum(!is.na(Alpha.df$inverse_simpson)))
#y<-as.data.frame(y)

#####to put on stats
summary(lm(inverse_simpson ~ as.factor(ACT) +Sex + age + batch , data = Alpha.df))


p<- ggboxplot(Alpha.df, x = "ACT", y = "inverse_simpson",##change depending on set of genes
              color = "ACT",
              add = "jitter",
              legend="")



p<-p+ ylab("Diversity") + xlab("") + ggtitle("Inverse Simpsons' diversity Index UC")+
  rotate_x_text(angle = 45)


#####add linear model pvals instead wilcox
stat.test <- compare_means(
  inverse_simpson ~ ACT, data = Alpha.df,
  method = "wilcox.test"
)

stat.test$p<-c(0.787,0.875,0.702)

stat.test <- stat.test %>%
  mutate(y.position = c(64, 60, 67))

p + stat_pvalue_manual(stat.test, label = "p")


write.csv(Alpha.df[,c("ACT","GSL.ID","inverse_simpson")],"Diversity_UC.csv")

###Differential Species

###DESeq approach
library(DESeq2)


#Run differential metafeature abundance analysis using DESeq2.

res.df<-res.df[,colnames(res.df) %in% samplesheet$GSL.ID]
samplesheet<-samplesheet[samplesheet$GSL.ID %in% colnames(res.df),]

#2. Estimate Dispersions
#ddsLovedis<-estimateDispersions(ddsLove)
#estimatedispersions<-as.data.frame(counts(ddsLovedis, normalized=TRUE))
geoMeans<-apply(res.df, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))

#RNAcount <- DESeqDataSetFromMatrix(countData = res.df, colData = samplesheet[,c("treat","Sex","age","batch")], 
#                                   design = ~  treat + Sex + age+batch)

###for dysbiosis, imid vs ctrl

samplesheet$imid<-droplevels(samplesheet$imid)
RNAcount <- DESeqDataSetFromMatrix(countData = res.df, colData = samplesheet[,c("imid","Sex","age","batch")], 
                                   design = ~  imid + Sex + age+batch)
ddsLove<-estimateSizeFactors(RNAcount, geoMeans=geoMeans)

res <- results(DESeq(ddsLove,fitType="local"), contrast = c('imid', 'PSA', 'CTRL'))


####HI vs CTRL
res<-as.data.frame(res)
res.sig<-res[res$pvalue<.1,]


####merge
res.sig<-as.data.frame(res.sig)
res.sig$V1<-rownames(res.sig)
res.sig<-merge(res.sig,species[,c("V1","ID")],by="V1")


write.csv(res.sig,"Res_PSA_CTRL.csv")
######

res$V1<-rownames(res)
res<-merge(res,species[,c("V1","ID")],by="V1")

res$SIG<-ifelse(res$pvalue<.05,"YES","NO")
vp<-ggplot(data=res,aes(x=log2FoldChange,y=-log10(pvalue)))+
  geom_point(aes(size=-log(pvalue),fill=SIG), colour="black",shape=21)+
  geom_vline(xintercept=0,linetype="dashed",col="red")+
  #geom_vline(xintercept=c(-0.5,0.5),linetype="dashed",col="blue")+
  #ggtitle("Ulcerative Colitis HI vs LOW") +
  theme(legend.position="none")+
  theme_minimal()


library(ggrepel)

vp<-vp+geom_text_repel(data=dplyr::filter(res, log2FoldChange < -0.4 & pvalue <.05 | log2FoldChange > 0.1 & pvalue <.05 ), aes(label=ID),point.padding = 2)

png("Volcano_PSA_CTRL.png",width = 1200, height = 880,)
vp+ ggtitle("Microbiota in Blood PSA vs CTRL") 
dev.off()

#Dysbiosis


####ad species to res.df
res.df$V1<-rownames(res.df)
res.df<-merge(res.df,species[,c("V1","ID")],by="V1")


#####species decrease 

mdi_dec <- unique(res.sig[res.sig$log2FoldChange<0,"ID"])
####avoid E.coli in down

mdi_dec<-mdi_dec[!(mdi_dec %in% "Escherichia_coli")]

mdi_down <- res.df[which(res.df$ID %in% mdi_dec), ]

#####species increase 

mdi_inc <- unique(res.sig[res.sig$log2FoldChange>0,"ID"])

mdi_up <- res.df[which(res.df$ID %in% mdi_inc), ]

A<-colSums(mdi_up[,!(colnames(mdi_up) %in% c("V1","ID"))]) + 1
B<-colSums(mdi_down[,!(colnames(mdi_down) %in% c("V1","ID"))]) + 1


Dysbiosis <- log(A/B)

Dysbiosis<-as.data.frame(Dysbiosis)
Dysbiosis$GSL.ID<-rownames(Dysbiosis)


library(ggsignif)


Alpha.df<-merge(Alpha.df,Dysbiosis,by="GSL.ID")



p<- ggboxplot(Alpha.df, x = "ACT", y = "Dysbiosis.y",##change depending on set of genes
              color = "ACT",
              add = "jitter",
              legend="")



p<-p+ ylab("Dysbiosis") + xlab("") + ggtitle("Dysbiosis PSA")+
  rotate_x_text(angle = 45)

#####add linear model pvals instead wilcox
stat.test <- compare_means(
  Dysbiosis ~ ACT, data = Alpha.df,
  method = "wilcox.test"
)


###Dysbiosis adjusted
summary(lm(Dysbiosis.y ~ as.factor(ACT) + Sex +age+ batch , data = Alpha.df))
###pvals coming from linear model

stat.test$p<-c(0.786,0.0397, 0.1146)


stat.test <- stat.test %>%
  mutate(y.position = c(12.7, 11, 12))
p + stat_pvalue_manual(stat.test, label = "p")


p<- p + stat_compare_means(method="wilcox.test",comparisons = my_comparisons,
                           #
                           #label = "p.signif",###only signifficance symbol stars
                           # method.args = list(var.equal = TRUE),##make t-test simila to ANOVA 2 groups,
                           
)
p


#The microbial dysbiosis index
#was calculated in R for each sample, according to Shaw et al. [5].
#Briefly, the microbial dysbiosis index
#was defined as log10 of the total abundance of taxa increased in CD divided 
#by taxa reduced in CD.
write.csv(Alpha.df[,c("ACT","GSL.ID","Dysbiosis")],"Disbiosis_PSA_based_res.csv")




