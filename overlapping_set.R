setwd("/mnt/md127/Analysis_Projects/CelGene_2018/")
load("CellGene_RNAseq_Counts_Salmon.RData")

sampleTable<-read.csv("Project4_IBD.RNAseq_Mastersheet.csv")
rownames(sampleTable)<-sampleTable$rna.code
sampleTable$cond<-paste(sampleTable$imid,sampleTable$activity,sep="_")

###select only samples with high qualiy, i.e high levels of RNA Integrity Number RIN >7

sampleTable<-sampleTable[sampleTable$rna.RIN > 7,]
sampleTable<-sampleTable[!(is.na(sampleTable$rna.RIN)),]
#Filtered out samples low corelation  :  IXTCB00176 IXTCB00031 IXTCB00178 IXTCB00244 IXTCB00388

##load interest genes

my.data<-read.delim("../Sanofi_2018/gwas_catalog_v1.0-associations_e91_r2018-01-16.csv",head=T)

#Per CD, usant la columna "Disease Trait" (uppercase): 
# ## [1] "CROHN'S DISEASE"                           
### [2] "CROHN'S DISEASE AND SARCOIDOSIS (COMBINED)"
## [3] "CROHN'S DISEASE AND PSORIASIS"             
## [4] "CROHN'S DISEASE AND CELIAC DISEASE"        
## [5] "ULCERATIVE COLITIS OR CROHN'S DISEASE"

#Per UC:
## [1] "ULCERATIVE COLITIS"                                      
## [2] "ULCERATIVE COLITIS OR CROHN'S DISEASE"                   
## [3] "SCLEROSING CHOLANGITIS AND ULCERATIVE COLITIS (COMBINED)"

#Agafar només loci a on "P.VALUE" < 5e-8.


####uc 

uc.genes<-my.data[grep("colitis",my.data$DISEASE.TRAIT),]
uc.genes$DISEASE.TRAIT<-droplevels(uc.genes$DISEASE.TRAIT)
uc.genes<-uc.genes[uc.genes$P.VALUE<5e-8,]
uc.genes<-uc.genes$REPORTED.GENE.S.
uc.genes<-unique(unlist(strsplit(as.character(uc.genes),",")))
uc.genes<-gsub(" ","",uc.genes)
extragenes.uc<-c("PRRCA","CXCR2","CXCR1","FA123B","PGAP3","CCR2","PUS10","UQCC","MIEN1")
#uc.genes[!(uc.genes %in% rownames(cts))]

uc.genes<-c(uc.genes,extragenes.uc)


###cd 

uc.genes<-my.data[grep("^Crohn's disease$|Ulcerative colitis or Crohn's disease|Crohn's disease and sarcoidosis|Crohn's disease and psoriasis|Crohn's disease and celiac disease",my.data$DISEASE.TRAIT),]
uc.genes$DISEASE.TRAIT<-droplevels(uc.genes$DISEASE.TRAIT)
uc.genes<-uc.genes[uc.genes$P.VALUE<5e-8,]
uc.genes<-uc.genes$REPORTED.GENE.S.
uc.genes<-unique(unlist(strsplit(as.character(uc.genes),",")))
uc.genes<-gsub(" ","",uc.genes)
extragenes.uc<-c("PRRCA","CXCR2","CXCR1","FA123B","PGAP3","CCR2","PUS10","UQCC","MIEN1")
#uc.genes[!(uc.genes %in% rownames(cts))]

uc.genes<-c(uc.genes,extragenes.uc)
# Following character vectors contain names of: 
# all.genes - all genes in dataset (17119) 
# de.genes - differentially expressed genes in dataset (86) 
# scs.genes - stem cell signature genes present in dataset (454)
# de.vs.scs - DE genes present in stem cell signature (9) 


## Run a Fisher's exact test on DE genes in SCS
# Generate a contingency table 

#pop size : Genome size, CD 12516; UC 12489
#sample size : CD 485; UC 383
#Number of items in the pop that are classified as successes : CD 1417; UC 2525
#Number of items in the sample that are classified as successes : CD 185;UC 146

#To compute a hypergeometric test, is that correct :
n.all.genes<-12516
n.de.genes<-1417
n.scs.genes<-485
n.de.vs.scs<-185

c.table <- matrix(c(n.all.genes, n.de.genes, n.scs.genes, n.de.vs.scs), nrow = 2) 

# Run Fisher's test and assign observed p value to 'z.obs' 
z.obs <- fisher.test(c.table)$p.value

## Write a function to generate a random sample of DE genes from all 
## the genes in array and calculate Fisher's statistics on it 

boot.f.test <- function(k){ # Where k is randomly drawn sample # of 86 genes from 17119 gene pool 
  
  n <- length(intersect(k, uc.genes)) # Compare sampled DE genes to the ones in Risk genes
  
  c.table <- matrix(c(n.all.genes, n.de.genes, n.scs.genes, n), nrow = 2) 
  
  fisher.test(c.table) }

## Resample DE genes 1000 times and calculate Fisher's exact test 

all.gene.names<-read.csv("../Sanofi_2018/Genes_all_CD_Sanofi.csv")
all.gene.names<-as.character(all.gene.names$x)
de.genes.names<-read.csv("filtered/CTRL_CD_filtered.csv")
de.genes.names<-as.character(de.genes.names$X)

z <- vector() 

for (i in 1:10000){ # Sample 86 genes from the pool of 17119 
  
  boot.de.genes <- sample(all.gene.names, size=length(de.genes.names), replace = T ) 
  
  # Calculate Fisher's test on the sample 
  
  z[i] <- boot.f.test(boot.de.genes)$p.value } 

# Calculate p value as a proportion of simulated p values equal or # smaller than observed p value

boot.p.value <- length(which(z <= z.obs))/length(z) 

boot.p.value

### End of Script

h<-hist(log(z), breaks=20,main=paste("Crohn's Disease",sep=" ","\np-value:","<0.0001"), xlab="log(pval)",col="lightblue",xlim=c(0,-80))

xfit<-seq(min(log(z)),max(log(z)),length=40)
yfit<-dnorm(xfit,mean=mean(log(z)),sd=sd(log(z)))
yfit <- yfit*diff(h$mids[1:2])*length(log(z))

lines(xfit, yfit, col="black", type="l",lwd=2,lty=2) 

abline(v=log(z.obs), lwd=2, col='red')
abline(v=mean(log(z)), lwd=2, col='black')
arrows(log(z.obs), 800, mean(log(z)), 800, col= 'darkgrey',lty=2,length=0.2)
arrows(mean(log(z)), 800, log(z.obs), 800, col= 'darkgrey',lty=2,length=0.2)
rug(log(z), col="blue")

library(Rmisc)

CI(log(z),ci=0.95)  





###aproach with uc.genes data only

exprs.uc.genes<-cpm_log[rownames(cpm_log) %in% uc.genes,]
exprs.uc.genes<-as.data.frame(colMeans(exprs.uc.genes))
exprs.uc.genes<-merge(exprs.uc.genes,samples2test[,c("imid","cond")],by="row.names")
exprs.uc.genes<-exprs.uc.genes[,c(2,3)]
colnames(exprs.uc.genes)<-c("Exprs","Group")

perm.t.test <- function (data,format,B=1000){
  options(scipen=999)
  if (format=="long") {
    unstacked.data <- unstack(data) #requires 'plyr'
    sample1 <- unstacked.data[[1]]
    sample2 <- unstacked.data[[2]]
  } else {
    sample1 <- data[,1]
    sample2 <- data[,2]
  }
  #get some statistics for the two samples
  n1 <- length(sample1)
  n2 <- length(sample2)
  mean1 <- round(mean(sample1), 2)
  mean2 <- round(mean(sample2),2)
  error1 <- qnorm(0.975)*sd(sample1)/sqrt(n1)
  error2 <- qnorm(0.975)*sd(sample2)/sqrt(n2)
  sample1_lci <- round(mean1 - error1,2)
  sample1_uci <- round(mean1 + error1,2)
  sample2_lci <- round(mean2 - error2,2)
  sample2_uci <- round(mean2 + error2,2)
  #get regular t-test results (equal variance)
  p.equal.var <- round(t.test(sample1, sample2, var.equal=TRUE)$p.value, 4)
  #get regular t-test results (unequal variance)
  p.unequal.var <- round(t.test(sample1, sample2, var.equal=FALSE)$p.value, 4)
  #start permutation procedures
  pooledData <- c(sample1, sample2)
  size.sample1 <- length(sample1)
  size.sample2 <- length(sample2)
  size.pooled <- size.sample1+size.sample2
  nIter <- B
  meanDiff <- numeric(nIter+1)
  meanDiff[1] <- round(mean1 - mean2, digits=2)
  for(i in 2:length(meanDiff)){
    index <- sample(1:size.pooled, size=size.sample1, replace=F)
    sample1.perm <- pooledData[index] 
    sample2.perm <- pooledData[-index] 
    meanDiff[i] <- mean(sample1.perm) - mean(sample2.perm) 
  }
  p.value <- round(mean(abs(meanDiff) >= abs(meanDiff[1])), digits=4)
  plot(density(meanDiff), main="Distribution of permuted mean differences", xlab="", sub=paste("sample 1 (n:", n1,") (95% CI lower bound., mean, 95% CI upper bound.):", sample1_lci, ",", mean1, ",", sample1_uci, "\nsample 2 (n:", n2,") (95% CI lower bound., mean, 95% CI upper bound.):", sample2_lci, ",", mean2, ",", sample2_uci,"\nobserved mean difference (dashed line):", meanDiff[1],"; permuted p.value (2-sided):", p.value, "( number of permutations:", B,")\nregular t-test p-values (2-sided):", p.equal.var, "(equal variance) ;", p.unequal.var, "(unequal variance)"), cex.sub=0.78)
  polygon(density(meanDiff), col="grey")
  rug(meanDiff, col="blue")
  abline(v=meanDiff[1], lty=2, col="red")
}


perm.t.test(exprs.uc.genes, format="long", B=10000)