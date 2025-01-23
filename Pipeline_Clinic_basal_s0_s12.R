library(tximport)
library(sva)
library(RUVSeq)
library(splines)
library(limma)
setwd("/Users/agomez/IMIDOMICS/PACTABA_Clinic/")

sampleTable<-read.csv("../PACTABA/PACTABA_RNAseq_Mastersheet.csv")
rownames(sampleTable)<-sampleTable$GSL.ID

codes<-read.csv("../PACTABA/Correspondencia_codis.csv")
sampleTable<-merge(sampleTable,codes,by="rna.code")
###pheno

pheno<-read.csv("../PACTABA/Pheno_sex_age.csv")
colnames(pheno)[3]<-"Code"
##2 is female
sampleTable<-merge(sampleTable,pheno,by="Code",all.x=T)
rownames(sampleTable)<-sampleTable$GSL.ID

week0<-read.csv("../PACTABA/Week0_PACTABA.csv")
week12<-read.csv("../PACTABA/Week12_PACTABA.csv")
week24<-read.csv("../PACTABA/Week24_PACTABA.csv")
week48<-read.csv("../PACTABA/Week48_PACTABA_F.csv")

####add 0 and 24 week
sampleTable<-merge(sampleTable,week0,by="Code",all.x=T)
#sampleTable<-merge(sampleTable,week24,by="Code",all.x=T)
sampleTable<-merge(sampleTable,week12,by="Code",all.x=T)


rownames(sampleTable)<-sampleTable$GSL.ID

library(ggpubr)
##########establish delta EULAR
sampleTable$delta<-sampleTable$das28esr_w0-sampleTable$das28esr_w12
sampleTable$EULAR<-ifelse(sampleTable$das28esr_w12<=3.2 & sampleTable$delta>1.2,"GOOD",
                          ifelse(sampleTable$das28esr_w12<=3.2 & sampleTable$delta>0.6 & sampleTable$delta<=1.2 ,"MOD",
                                 ifelse(sampleTable$das28esr_w12<=3.2 & sampleTable$delta<=0.6,"NON_RESP",
                                        ifelse(sampleTable$das28esr_w12>3.2 & sampleTable$das28esr_w12<=5.1 & sampleTable$delta>1.2,"MOD",
                                               ifelse(sampleTable$das28esr_w12>3.2 &  sampleTable$das28esr_w12<=5.1 & sampleTable$delta>0.6 & sampleTable$delta<=1.2,"MOD",
                                                      ifelse(sampleTable$das28esr_w12>3.2 & sampleTable$das28esr_w12<=5.1 & sampleTable$delta<=0.6,"NON_RESP",
                                                             ifelse(sampleTable$das28esr_w12>5.1 & sampleTable$delta>1.2,"MOD",
                                                                    ifelse(sampleTable$das28esr_w12> 5.1 & sampleTable$delta>0.6 & sampleTable$delta<=1.2 ,"NON_RESP",
                                                                           ifelse(sampleTable$das28esr_w12> 5.1 & sampleTable$delta<=0.6,"NON_RESP","NA")))))))))###############EULAR

##write out pheno PACTABA to check
pheno_PACTABA<-sampleTable[,c("Code","das28esr_w0","das28esr_w12","delta","EULAR","donor.code.x")]

pheno_PACTABA<-pheno_PACTABA[!(is.na(pheno_PACTABA$delta)),]

##change iD
pheno_PACTABA$Plasmes<-gsub("PACT_0","",pheno_PACTABA$donor.code.x)
pheno_PACTABA$Plasmes<-gsub("^0","",pheno_PACTABA$Plasmes)
pheno_PACTABA$Plasmes<-paste(pheno_PACTABA$Plasmes,"-S0",sep="")###change S0, S12, S24 depending on time analyzed
pheno_PACTABA<-pheno_PACTABA[,c("Plasmes","delta","EULAR")]
pheno_PACTABA<-unique(pheno_PACTABA,MARGIN=1)

sampleTable$Plasmes<-gsub("PACT_0","",sampleTable$donor.code.x)
sampleTable$Plasmes<-gsub("^0","",sampleTable$Plasmes)
sampleTable$Plasmes<-paste(sampleTable$Plasmes,"-S0",sep="")
acCARP<-merge(acCARP,sampleTable,by="Plasmes")
acCARP$wk<-ifelse(grepl("-S0",acCARP$donationId),"wk0",ifelse(grepl("-S12",acCARP$donationId),"wk12",
                   ifelse(grepl("-S24",acCARP$donationId),"wk24","wk48"                         
 
                                                                                     )))


sel<-acCARP[acCARP$wk=="wk0"|acCARP$wk=="wk12", ]


my.samples<-c( "5697-RMM-0002", "5697-RMM-0001", "5697-RMM-0008", "5697-RMM-0007", "5697-RMM-0045", "5697-RMM-0046", "5697-RMM-0112", "5697-RMM-0113",
               "5697-RMM-0140" ,"5697-RMM-0141" ,"5697-RMM-0168", "5697-RMM-0167", "5697-RMM-0189", "5697-RMM-0188", "5697-RMM-0135", "5697-RMM-0134",
               "5697-RMM-0178", "5697-RMM-0177" ,"5697-RMM-0153", "5697-RMM-0152", "5697-RMM-0198", "5697-RMM-0197", "5697-RMM-0128", "5697-RMM-0129",
               "5697-RMM-0031", "5697-RMM-0030" ,"5697-RMM-0097", "5697-RMM-0096", "5697-RMM-0015", "5697-RMM-0014", "5697-RMM-0025", "5697-RMM-0026",
               "5697-RMM-0101" ,"5697-RMM-0102" ,"5697-RMM-0056", "5697-RMM-0055", "5697-RMM-0108", "5697-RMM-0107", "5697-RMM-0087", "5697-RMM-0086",
               "5697-RMM-0066" ,"5697-RMM-0067" ,"5697-RMM-0035", "5697-RMM-0036", "5697-RMM-0062", "5697-RMM-0061", "5697-RMM-0020", "5697-RMM-0021",
               "5697-RMM-0071" ,"5697-RMM-0070" ,"5697-RMM-0076", "5697-RMM-0075", "5697-RMM-0080", "5697-RMM-0081", "5697-RMM-0123", "5697-RMM-0124",
               "5697-RMM-0158" ,"5697-RMM-0157" ,"5697-RMM-0204", "5697-RMM-0203", "5697-RMM-0040", "5697-RMM-0041", "5697-RMM-0183", "5697-RMM-0182",
               "5697-RMM-0194" ,"5697-RMM-0193" ,"5697-RMM-0091", "5697-RMM-0090", "5697-RMM-0119", "5697-RMM-0118", "5697-RMM-0173", "5697-RMM-0172",
               "5697-RMM-0163" ,"5697-RMM-0164" ,"5697-RMM-0146", "5697-RMM-0147")

sel<-sel[sel$GSL.ID %in% my.samples,]

to.table<-sel[sel$wk=="wk0",]
to.table<-sel[sel$wk=="wk12",]


pheno2<-read.delim("Extra_pheno_info.csv")
to.table2<-pheno2[pheno2$ID %in% to.table$Code,]
###get only wk0 and wk12

###Threshold acCARP 132.5
acCARP<-read.csv("Plasmes_S0.csv")
acCARP<-read.csv("Plasmes_S12.csv")

####in case acCARP at w12, change S0 for S12
#pheno_PACTABA$Plasmes<-gsub("S0","S12",pheno_PACTABA$Plasmes)

acCARP<-merge(acCARP,pheno_PACTABA,by="Plasmes")###41 samples w24, 57 w12
acCARP$CARP<-ifelse(acCARP$UA.mL>132.5,"POS","NEG")
acCARP$RESP<-ifelse(acCARP$EULAR=="NON_RESP","NON_RESP","RESP")

###delta DAS Plot
sampleTable<-read.csv("../PACTABA/PACTABA_RNAseq_Mastersheet.csv")
rownames(sampleTable)<-sampleTable$GSL.ID

codes<-read.csv("../PACTABA/Correspondencia_codis.csv")
sampleTable<-merge(sampleTable,codes,by="rna.code")
###pheno

pheno<-read.csv("../PACTABA/Pheno_sex_age.csv")
colnames(pheno)[3]<-"Code"
##2 is female
sampleTable<-merge(sampleTable,pheno,by="Code",all.x=T)
rownames(sampleTable)<-sampleTable$GSL.ID

sampleTable<-merge(sampleTable,week0,by="Code",all.x=T)
sampleTable<-merge(sampleTable,week12,by="Code",all.x=T)

###to.plot

to.plot<-sampleTable[,c("donor.code.x","das28esr_w0","das28esr_w12")]
to.plot<-to.plot[complete.cases(to.plot),]

to.plot$Plasmes<-gsub("PACT_0","",to.plot$donor.code.x)
to.plot$Plasmes<-gsub("^0","",to.plot$Plasmes)
to.plot$Plasmes<-paste(to.plot$Plasmes,"-S0",sep="")###change S0, S12, S24 depending on time analyzed
acCARP<-merge(acCARP,to.plot,by="Plasmes")###41 samples
colnames(acCARP)<-gsub(".x","",colnames(acCARP))
acCARP$RESP<-ifelse(acCARP$EULAR=="NON_RESP","NON_RESP","RESP")


table(acCARP$RESP) #23 NON 34 YES 
##put in another shape
##test
##DAS28_Timecourse_CARP

library(reshape)
mdata <- melt(acCARP[, c("RESP","das28esr_w0","das28esr_w12")], id=c("RESP"))

###subset POS responders and NEG non responders
to.test<-acCARP[(acCARP$CARP=="POS" & acCARP$RESP=="RESP") | (acCARP$CARP=="NEG" & acCARP$RESP=="NON_RESP") ,]## 18 resp Non 12
to.test<-acCARP[(acCARP$ACPA=="POS" & acCARP$RESP=="RESP") | (acCARP$ACPA=="NEG" & acCARP$RESP=="RESP") ,]## 18 resp Non 12

mdata <- melt(acCARP[, c("CARP","das28esr_w0","das28esr_w12")], id=c("CARP"))
mdata <- melt(acCARP[, c("ACPA","das28esr_w0","das28esr_w12")], id=c("ACPA"))
mdata <- melt(to.test[, c("ACPA","das28esr_w0","das28esr_w12")], id=c("ACPA"))###Fig 1
mdata <- melt(to.test[, c("CARP","das28esr_w0","das28esr_w12")], id=c("CARP"))

mdata$variable<-gsub("das28esr\\_","",mdata$variable)

colnames(mdata)[1]<-"anti_CarP"
#https://stats.stackexchange.com/questions/126375/how-to-determine-signficant-difference-between-2-curves
sum<-lm(value~variable+RESP+variable*RESP, data=mdata)
sum<-lm(value~variable+ACPA+variable*ACPA, data=mdata)

##get pval

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

pval<-lmp(sum)
pval<-signif(pval, digits=3)

tiff("Fig1c.tiff", units="in", width=5, height=5, res=750)

p<-ggline(mdata, x = "variable", y = "value",color = "anti_CarP",add = c("mean_se"))
p<-ggline(mdata, x = "variable", y = "value",color = "ACPA",add = c("mean_se"))

p<-p+ ylab("DAS28") + xlab("") + ggtitle(paste("pval=",pval,sep=""))

p

dev.off()

###3rd question
### Remission patients with high CARP concentration?? 

my.data<-acCARP[,c("UA.mL","RESP")]
my.data<-acCARP[,c("CCPLUS","RESP")]

##s0
my.data<-acCARP[grepl("S0",acCARP$Plasmes,perl=T),]
my.data<-my.data[,c("UA.mL","RESP")]

##s12

my.data<-acCARP[grepl("S12",acCARP$Plasmes,perl=T),]
my.data<-my.data[,c("UA.mL","RESP")]

cols<-c("lightblue","red")


p <- ggboxplot(my.data, x = "RESP", y = "CCPLUS",##change depending on set of genes
               color = "RESP",
               add = "jitter",
               legend="")


p<-p+ ggtitle(paste("EULAR RESPONSE ACPA week 0"))+
  rotate_x_text(angle = 45)
# Specify the comparisons you want
my_comparisons <- list( c("RESP", "NON_RESP"))

p<- p + stat_compare_means(method="t.test",comparisons = my_comparisons,
                           #
                           #label = "p.signif",###only signifficance symbol stars
                           method.args = list(var.equal = TRUE),##make t-test simila to ANOVA 2 groups,
                           symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                               symbols = c("****", "***", "**", "*", "ns"))
)

######Response associated to decrease in anti_CARP []


######anti-CARP []

carp0<-read.csv("Plasmes_S0.csv")
carp12<-read.csv("Plasmes_S12.csv")

carp<-rbind(carp0,carp12)
carp$Time<-ifelse(grepl("S0",carp$Plasmes),"w0","w12")
carp$TIME<-ifelse(grepl("S0",carp$Plasmes),"0","12")
carp$TIME<-as.numeric(carp$TIME)
carp$ID<-gsub("S0|S12","PACT",carp$Plasmes)
##Response


sampleTable$RESP<-ifelse(sampleTable$EULAR=="NON_RESP","NON_RESP","RESP")
rems<-sampleTable[,c("donor.code.x","RESP")]
rems$ID<-gsub("PACT_0","",rems$donor.code.x)
rems$ID<-gsub("^0","",rems$ID)
rems$ID<-paste(rems$ID,"-PACT",sep="")
rems<-unique(rems,MARGIN=1)
rems<-rems[!(is.na(rems$RESP)),]

carp<-merge(carp,rems, by="ID")
carp<-carp[!(is.na(carp$RESP)),]
carp$Time<-as.factor(carp$Time)
###to.plot

colnames(carp)<-gsub(".x","",colnames(carp))

carp$CARP<-ifelse(carp$UA.mL>132.5,"POS","NEG")
##subset wk12

wk12<-carp[carp$Time=="w12",]
wk0<-carp[carp$Time=="w0",]

table(carp[,c("RESP","CARP")])

##3test of the slopes

sum<-lm(TIME~UA.mL*RESP, data=carp)
sum<-lm(UA.mL~TIME+RESP+TIME*RESP, data=carp)

####test of slopes

#fit <- lm(y ~ group + x + x:group) 
#where y is the response of the 2 groups.
#The p-value of x:group gives the probability for the two slopes to be   
#different, and the estimated values of parameters are these of both   
#populations. 

#+ggtitle()

####test of equality of slopes

sum<-lm(UA.mL~RESP + Time+ Time*RESP, data = carp)

anova(sum)###get resp


lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

pval<-lmp(sum)
pval<-signif(pval, digits=3)

pval<-"0.02373"
###Fig 1b

tiff("Fig1b.tiff", units="in", width=5, height=5, res=750)

p<-ggline(carp, x = "Time", y = "UA.mL",color = "RESP",add = c("mean_se"))+theme(legend.title=element_blank())

p<-p+ ylab("anti-CarP Level") + xlab("") + ggtitle(paste("pval=",pval))
p

dev.off()

##means


Resp<-acCARP[acCARP$RESP=="RESP",]
NResp<-acCARP[!(acCARP$RESP=="RESP"),]


Resp<-carp[carp$RESP=="RESP",]
NResp<-carp[!(carp$RESP=="RESP"),]

mean(Resp[Resp$Time=="w12","UA.mL"])





########################################SAME ANALYSIS with ACPA


sampleTable<-read.csv("../PACTABA/PACTABA_RNAseq_Mastersheet.csv")
rownames(sampleTable)<-sampleTable$GSL.ID

codes<-read.csv("../PACTABA/Correspondencia_codis.csv")
sampleTable<-merge(sampleTable,codes,by="rna.code")
###pheno

pheno<-read.csv("../PACTABA/Pheno_sex_age.csv")
colnames(pheno)[3]<-"Code"
##2 is female
sampleTable<-merge(sampleTable,pheno,by="Code",all.x=T)
rownames(sampleTable)<-sampleTable$GSL.ID

week0<-read.csv("../PACTABA/Week0_PACTABA.csv")
week12<-read.csv("../PACTABA/Week12_PACTABA.csv")
week24<-read.csv("../PACTABA/Week24_PACTABA.csv")
week48<-read.csv("../PACTABA/Week48_PACTABA_F.csv")
####add 0 and 24 week
sampleTable<-merge(sampleTable,week0,by="Code",all.x=T)
#sampleTable<-merge(sampleTable,week24,by="Code",all.x=T)
sampleTable<-merge(sampleTable,week12,by="Code",all.x=T)


rownames(sampleTable)<-sampleTable$GSL.ID

##########establish delta EULAR
sampleTable$delta<-sampleTable$das28esr_w0-sampleTable$das28esr_w12
sampleTable$EULAR<-ifelse(sampleTable$das28esr_w12<=3.2 & sampleTable$delta>1.2,"GOOD",
                          ifelse(sampleTable$das28esr_w12<=3.2 & sampleTable$delta>0.6 & sampleTable$delta<=1.2 ,"MOD",
                                 ifelse(sampleTable$das28esr_w12<=3.2 & sampleTable$delta<=0.6,"NON_RESP",
                                        ifelse(sampleTable$das28esr_w12>3.2 & sampleTable$das28esr_w12<=5.1 & sampleTable$delta>1.2,"MOD",
                                               ifelse(sampleTable$das28esr_w12>3.2 &  sampleTable$das28esr_w12<=5.1 & sampleTable$delta>0.6 & sampleTable$delta<=1.2,"MOD",
                                                      ifelse(sampleTable$das28esr_w12>3.2 & sampleTable$das28esr_w12<=5.1 & sampleTable$delta<=0.6,"NON_RESP",
                                                             ifelse(sampleTable$das28esr_w12>5.1 & sampleTable$delta>1.2,"MOD",
                                                                    ifelse(sampleTable$das28esr_w12> 5.1 & sampleTable$delta>0.6 & sampleTable$delta<=1.2 ,"NON_RESP",
                                                                           ifelse(sampleTable$das28esr_w12> 5.1 & sampleTable$delta<=0.6,"NON_RESP","NA")))))))))###############EULAR

##write out pheno PACTABA to check
pheno_PACTABA<-sampleTable[,c("Code","das28esr_w0","das28esr_w12","delta","EULAR","donor.code.x")]

pheno_PACTABA<-pheno_PACTABA[!(is.na(pheno_PACTABA$delta)),]


##########establish delta EULAR


##write out pheno PACTABA to check
pheno_PACTABA<-sampleTable[,c("Code","rna.code","donor.code.x","donationId","das28esr_w0","das28esr_w12","delta","EULAR")]

pheno_PACTABA<-pheno_PACTABA[!(is.na(pheno_PACTABA$delta)),]

##change iD
pheno_PACTABA$Plasmes<-gsub("PACT_0","",pheno_PACTABA$donor.code.x)
pheno_PACTABA$Plasmes<-gsub("^0","",pheno_PACTABA$Plasmes)
pheno_PACTABA$Plasmes<-paste(pheno_PACTABA$Plasmes,"-S0",sep="")###change S0, S12, S24 depending on time analyzed
pheno_PACTABA<-pheno_PACTABA[,c("Plasmes","delta","EULAR")]
pheno_PACTABA<-unique(pheno_PACTABA,MARGIN=1)
###Threshold acCARP 132.5
acCARP<-read.csv("ACPA_S0.csv")
acCARP<-merge(acCARP,pheno_PACTABA,by="Plasmes")###41 samples w24, 57 w12
acCARP$ACPA<-ifelse(acCARP$CCPLUS>50,"POS","NEG")
acCARP$RESP<-ifelse(acCARP$EULAR=="NON_RESP","NON_RESP","RESP")

########

#Boxplot deltas

###
my.data<-acCARP[,c("delta","ACPA")]
cols<-c("red","lightblue")

tiff("Fig1c_Boxplot_deltas.tiff", units="in", width=5, height=5, res=750)


p <- ggboxplot(my.data, x = "ACPA", y = "delta",##change depending on set of genes
               color = "ACPA",
               add = "jitter",
               legend="")


p<-p +rotate_x_text(angle = 45)+ylab("deltaDAS")+ theme(axis.title.x = element_blank()) + ggtitle(paste("pval=0.45"))

# Specify the comparisons you want
my_comparisons <- list( c("POS", "NEG"))

p<- p + stat_compare_means(method="t.test",comparisons = my_comparisons,
                           #
                           #label = "p.signif",###only signifficance symbol stars
                           #method.args = list(var.equal = TRUE),##make t-test simila to ANOVA 2 groups,
                           symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                               symbols = c("****", "***", "**", "*", "ns"))
)

p

dev.off()



#######################
##Fig 1d
#my.data<-acCARP[,c("delta","ACPA")]
my.data<-acCARP[,c("CCPLUS","RESP")]
colnames(my.data)[1]<-"ACPA"
cols<-c("red","lightblue")

tiff("Fig1d.tiff", units="in", width=5, height=5, res=750)


p <- ggboxplot(my.data, x = "RESP", y = "ACPA",##change depending on set of genes
               color = "RESP",
               add = "jitter",
               legend="")


p<-p +rotate_x_text(angle = 45)+ylab("ACPA Level")+ theme(axis.title.x = element_blank()) + ggtitle(paste("pval=0.28"))
 
# Specify the comparisons you want
my_comparisons <- list( c("NON_RESP", "RESP"))

p<- p + stat_compare_means(method="t.test",comparisons = my_comparisons,
                           #
                           #label = "p.signif",###only signifficance symbol stars
                           #method.args = list(var.equal = TRUE),##make t-test simila to ANOVA 2 groups,
                           symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                               symbols = c("****", "***", "**", "*", "ns"))
)

p

dev.off()


###Fig 1d acCARP instead ACPA, same procedure as Fig 1d

###Threshold acCARP 132.5
acCARP<-read.csv("Plasmes_S0.csv")
acCARP<-merge(acCARP,pheno_PACTABA,by="Plasmes")###41 samples w24, 57 w12
#acCARP$ACPA<-ifelse(acCARP$CCPLUS>50,"POS","NEG")
acCARP$CARP<-ifelse(acCARP$UA.mL>132.5,"POS","NEG")

acCARP$RESP<-ifelse(acCARP$EULAR=="NON_RESP","NON_RESP","RESP")

#######################
##Fig 1d
#my.data<-acCARP[,c("delta","ACPA")]
my.data<-acCARP[,c("UA.mL","RESP")]
colnames(my.data)[1]<-"acCARP"
cols<-c("red","lightblue")

tiff("Fig1b_nueva_acCARP_basal_week.tiff", units="in", width=5, height=5, res=750)


p <- ggboxplot(my.data, x = "RESP", y = "acCARP",##change depending on set of genes
               color = "RESP",
               add = "jitter",
               legend="")


p<-p +rotate_x_text(angle = 45)+ylab("acCARP Level")+ theme(axis.title.x = element_blank()) + ggtitle(paste("pval=0.018"))

# Specify the comparisons you want
my_comparisons <- list( c("NON_RESP", "RESP"))

p<- p + stat_compare_means(method="t.test",comparisons = my_comparisons,
                           #
                           #label = "p.signif",###only signifficance symbol stars
                           #method.args = list(var.equal = TRUE),##make t-test simila to ANOVA 2 groups,
                           symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                               symbols = c("****", "***", "**", "*", "ns"))
)

p

dev.off()




###question 1

##delta DAS is different between CARP + and CARP -
##delta DAS is different between ACPA + and ACPA -



###delta DAS Plot
sampleTable<-read.csv("../PACTABA/PACTABA_RNAseq_Mastersheet.csv")
rownames(sampleTable)<-sampleTable$GSL.ID

codes<-read.csv("../PACTABA/Correspondencia_codis.csv")
sampleTable<-merge(sampleTable,codes,by="rna.code")
###pheno

pheno<-read.csv("../PACTABA/Pheno_sex_age.csv")
colnames(pheno)[3]<-"Code"
##2 is female
sampleTable<-merge(sampleTable,pheno,by="Code",all.x=T)
rownames(sampleTable)<-sampleTable$GSL.ID

sampleTable<-merge(sampleTable,week0,by="Code",all.x=T)
sampleTable<-merge(sampleTable,week12,by="Code",all.x=T)
#sampleTable<-merge(sampleTable,week24,by="Code",all.x=T)

###to.plot

to.plot<-sampleTable[,c("donor.code.x","das28esr_w0","das28esr_w12")]
to.plot<-unique(to.plot,MARGIN=1)

to.plot$Plasmes<-gsub("PACT_0","",to.plot$donor.code.x)
to.plot$Plasmes<-gsub("^0","",to.plot$Plasmes)
to.plot$Plasmes<-paste(to.plot$Plasmes,"-S0",sep="")###change S0, S12, S24 depending on time analyzed
acCARP<-merge(acCARP,to.plot,by="Plasmes")###41 samples
colnames(acCARP)<-gsub(".x","",colnames(acCARP))


##put in another shape
##test


library(reshape)
mdata <- melt(acCARP[, c("ACPA","das28esr_w0","das28esr_w12")], id=c("ACPA"))
mdata$variable<-gsub("das28esr\\_","",mdata$variable)
#https://stats.stackexchange.com/questions/126375/how-to-determine-signficant-difference-between-2-curves
sum<-lm(value~variable+ACPA+variable*ACPA, data=mdata)
##get pval

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

pval<-lmp(sum)
pval<-signif(pval, digits=3)

p<-ggline(mdata, x = "variable", y = "value",color = "ACPA",add = c("mean_se"))

p<-p+ ylab("DAS28") + xlab("") + ggtitle(paste("DAS28 differences across time for ACPA groups at w0 pval=",pval,sep=""))




###2nd question
###Remission

###categorize using remission at w24
acCARP$REM<-ifelse(acCARP$das28esr_w24<3.2,"YES","NO")
#acCARP$REM<-ifelse(acCARP$das28esr_w12<3.2,"YES","NO")


###ACPA POS REM POS

table(acCARP[acCARP$ACPA=="POS","REM"])
table(acCARP[acCARP$ACPA=="NEG","REM"])

#####    ACPA

###REM     POS   NEG
##   YES    7    3
##   NO    25     6


ACPA_REM <-
  matrix(c(6, 21, 3, 4),
         nrow = 2,
         dimnames = list(REM = c("YES", "NO"),
                         ACPA = c("POS", "NEG")))
fisher.test(ACPA_REM)


###p-value = 0.6622 w0
#p-value = 1 w12
#p-value 0.3482 w24


###3rd question
### Remission patients with high CARP concentration?? 

my.data<-acCARP[,c("CCPLUS","REM")]
cols<-c("red","lightblue")


p <- ggboxplot(acCARP, x = "REM", y = "CCPLUS",##change depending on set of genes
               color = "REM",
               add = "jitter",
               legend="")


p<-p+ ggtitle(paste("Low activity/Remission ACPA week 12"))+
  rotate_x_text(angle = 45)
# Specify the comparisons you want
my_comparisons <- list( c("YES", "NO"))

p<- p + stat_compare_means(method="t.test",comparisons = my_comparisons,
                           #
                           #label = "p.signif",###only signifficance symbol stars
                           method.args = list(var.equal = TRUE),##make t-test similar to ANOVA 2 groups,
                           symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                               symbols = c("****", "***", "**", "*", "ns"))
)

p
######Response associated to decrease in anti_CARP []


######anti-CARP []

acpa0<-read.csv("PACTABA_Clinic/ACPA_S0.csv")
acpa12<-read.csv("PACTABA_Clinic/ACPA_S12.csv")

acpa<-rbind(acpa0,acpa12)
acpa$Time<-ifelse(grepl("S0",acpa$Plasmes),"w0","w12")
acpa$TIME<-ifelse(grepl("S0",acpa$Plasmes),"0","12")
acpa$TIME<-as.numeric(acpa$TIME)
acpa$ID<-gsub("S0|S12","PACT",acpa$Plasmes)
##Response
sampleTable$RESP<-ifelse(sampleTable$EULAR=="NON_RESP","NON_RESP","RESP")
rems<-sampleTable[,c("donor.code.x","RESP")]
rems$ID<-gsub("PACT_0","",rems$donor.code.x)
rems$ID<-gsub("^0","",rems$ID)
rems$ID<-paste(rems$ID,"-PACT",sep="")
rems<-unique(rems,MARGIN=1)
rems<-rems[!(is.na(rems$RESP)),]

acpa<-merge(acpa,rems, by="ID")
acpa<-acpa[!(is.na(acpa$RESP)),]
acpa$Time<-as.factor(acpa$Time)
###to.plot

colnames(acpa)<-gsub(".x","",colnames(acpa))
##3test of the slopes

sum<-lm(CCPLUS~TIME*RESP, data=acpa)
sum<-lm(UA.mL~RESP + Time + RESP:Time, data = acpa)
anova(sum)
#anova(lm(UA.mL~TIME,data=carp),lm(UA.mL~TIME*REM,data=carp))
##get pval

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

pval<-lmp(sum)
pval<-signif(pval, digits=3)

p<-ggline(acpa, x = "Time", y = "CCPLUS",color = "RESP",add = c("mean_se"))

p<-p+ ylab("[ACPA]") + xlab("") + ggtitle(paste("ACPA Conc Evolution pval <0.05"))
p

######


Resp<-acpa[acpa$RESP=="RESP",]
NResp<-acpa[!(acpa$RESP=="RESP"),]

mean(Resp[Resp$Time=="w12","CCPLUS"])
mean(Resp[Resp$Time=="w0","CCPLUS"])


mean(NResp[NResp$Time=="w12","CCPLUS"])
mean(NResp[NResp$Time=="w0","CCPLUS"])



