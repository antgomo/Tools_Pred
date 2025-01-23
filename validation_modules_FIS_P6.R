setwd("/mnt/md127/Analysis_Projects/FIS_Meth/")
library(igraph)
#Where spinmods is a list, with each element of the list a vector of Entrez IDs corresponding to the nodes
#of a community, pin is an iGraph graph object with nodes corrresponding to genes labelled by Entrez IDs,
#stat is a vector of node weights to be randomly permuted over the network and nrand is the number of
#random permutations to compute. In order to validate the T1D communities in the SKIN data set, it can be
#called as follows:




#However, the validation of a FEM module that you infer in some other dataset should not be too difficult.

#The way we had implemented this, if I recall correctly, was by comparing the weighted average statistic of DM for the module to the null 
#obtained by selecting statistics at random from the study and allocating them randomly to the nodes of the module. 
#The weights in the average do not change and would be +1 or -1 depending on the directionality of the observed statistic. 
#You would run say 1000 Monte-Carlo randomizations to get the null, from which you can then assess significance. Hope this helps. 

##############################3
##FUNCTIONS
#grap weighted according its methylation stats

statLabel = function(g, stat)
  # Label graph g with statistics stat (averaged on edges)
{
  A = get.adjacency(g); 
  A = A[names(stat),names(stat)]; # just to make sure
  astat = abs(stat);
  TMAX = max(astat);
  temp1 = apply(A, 1, function(v) return(v*astat))
  W = (temp1 + t(temp1))/(2*TMAX);
  print("Generating weighted graph");
  G = graph.adjacency(W, mode = "undirected", weighted=TRUE);
  V(G)$weight = astat;
  return(G);
}

Modularity = function(v, g)
  # Returns modularity score of subgraph induced by vertices v in iGraph g
{ 
  h = induced.subgraph(g, v);
  return(sum( E(h)$weight ))
}

nodeValidate = function(lv, g, tstats, nrandomizations = 1000)
  # Test elements of lv in g over node-randomizations of tstats 
{
  #A = get.adjacency(g);
  nrm.lv = list(); 
  nrp.v = vector();
  obs.v = vector();
  output = list();
  TMAX = max(tstats);
  
  # First restrict to g:
  mods = lapply(lv, function(v) {
    return( intersect(v,V(g)$name) )
  });
  print(sapply(mods, length));
  
  # Observed modularity values
  obs.v = as.vector(lapply(mods, function (j) { return(Modularity(j, g)) }));
  
  # Node-randomizations of modularity
  nmods = length(mods);
  nrm.lv = lapply(1:nmods, function (i) {
    j = mods[[i]];
    v = vector();
    h = induced.subgraph(g,j); # new iGraph uses induced.subgraph() not subgraph()
    B = get.adjacency(h, sparse=FALSE); # Do not want sparse object
    for (k in 1:nrandomizations) {
      print(paste("Testing significance: module",i,"of",nmods,"Randomization",k,"of",nrandomizations))
      atperm = sample(tstats, nrow(B) , replace=FALSE)
      temp1 = apply(B, 1, function(v) return(v*atperm))
      W = (temp1 + t(temp1))/(2*TMAX);
      v[k] = sum(W)/2; # W is weighted adj matrix (with diag=0) so every edge counted twice
      #v[k]<-sum(W)
    }
    return(v);
  }); names(nrm.lv) = names(mods);
  
  # Empirical p-values
  for (j in 1:length(mods)) {
    nrp.v[j] = length( which(nrm.lv[[j]] > obs.v[j]) )/nrandomizations;
  }; names(nrp.v) = names(mods);
  
  output[[1]] = nrp.v; output[[2]] = obs.v; output[[3]] = nrm.lv; 
  names(output) = c("fdr","Observed","Random");
  return(output);
}

# Code to assign colors to vectors
vect2color = function(v, palettev, breaks) {
  w = v; 
  for (i in 1:length(palettev)) { w[which( v >= breaks[i] & v < breaks[i+1] )] = palettev[i]; }
  return(w);
}

renderModule = function(eid, g, pval, stat) 
  # Renders module with ENTREZ IDs eid, on top of network g with pvalues pval and statistics stat:
{
  # Vertex palette:
  vertexPalette.v = maPalette(low="yellow", high= "blue", mid="white", k=49);
  vertexBreaks.v = seq(from = log10(0.05), to = -log10(0.05), 
                       by = -2*log10(0.05)/(-2+length(vertexPalette.v )))
  vertexBreaks.v = c(-(1+max(-log10(pval))), vertexBreaks.v, 1+max(-log10(pval)));
  length(vertexPalette.v); length(vertexBreaks.v);
  # Edge palette: grey to red via pink, red at top 10% of edges
  edgePalette.v = maPalette(low="lightgrey",high="red",mid="pink",k = 50);
  netPercentiles.v = seq(from=0,to=1,by=1/length(edgePalette.v))
  edgeBreaks.v = sapply(netPercentiles.v, function(x) quantile(E(g)$weight, x))
  edgeBreaks.v[1] = 0; edgeBreaks.v[length(edgeBreaks.v)] = 1;
  
  # Compute iGraph object
  h = induced.subgraph(g, eid);
  stat.v = stat[V(h)$name];
  pval.v = pval[V(h)$name];
  slpval.v = sign(stat.v)*-log10(pval.v); # signed, logged p-values
  
  par(mar=c(4,0,2,0))
  # Color edges between grey and red according to significance  
  E(h)$color = vect2color(E(h)$weight, edgePalette.v, edgeBreaks.v);
  # Color nodes blue to yellow according to hyper/hypo-methylation
  V(h)$color = vect2color(slpval.v, vertexPalette.v, vertexBreaks.v);
  
  if ( length(V(h)) <= 20 ) {  vl = unlist(entrez2symbol[V(h)$name])  } else { vl = "" }
  
  plot(h,layout = layout.fruchterman.reingold,  # same layout each time
       vertex.label = vl, 
       vertex.frame.color = "black", 
       vertex.label.dist = 1.1, 
       vertex.label.font = 3, vertex.label.color = "black", 
       #vertex.size = 15*13/length(V(h)),
       #edge.width = 160/length(V(h))
       vertex.size = if (length(V(h)) < 50) { 15 } else {15*13/length(V(h))},
       edge.width  = if (length(V(h)) < 50) {  6 } else { 160/length(V(h) )}
  );
}

####################################################################ANALYSIS
##load set to test P6
  

load("Validations/Epimods_analysis/Set1/Quant_baseline/EpiMod_wk0_DAS28.RData")#Toci in
load("/mnt/md127/Analysis_Projects/FIS_Meth/Results_antiTNF/delta/Epimod_Delta_wk0.RData")


Spain.Epi<-EpiMod.o

load("/mnt/md127/Analysis_Projects/FIS_Meth/Results_antiTNF/delta/P6_val/Epimod_Delta_wk0.RData")
load("/mnt/md127/Analysis_Projects/FIS_Meth/Validations/Epimods_analysis/Set2/Quant_baseline/Epimod_DAS28_wk0.RData")


##now get module names to test


lv<-lapply(Spain.Epi$topmod, function(x) {
  
  
  return(as.character(x[1]$EntrezID))
} 

)

####P6's graph to allocate randomly the stats

g<-graph_from_adjacency_matrix(EpiMod.o$adj)


##methylation stats to allocate them randomly 
load("/mnt/md127/Analysis_Projects/FIS_Meth/Results_antiTNF/delta/P6_val/intEPI_stats_P6_Delta_wk0.RData")
load("/mnt/md127/Analysis_Projects/FIS_Meth/Validations/Epimods_analysis/Set2/Quant_baseline/int_stats_DAS28_wk0.RData")

tstats<-intEpi.o$statM$t
names(tstats)<-rownames(intEpi.o$statM)


##build graph

g<-statLabel(g, tstats)

##run algorithm

Stats<-nodeValidate(lv,g, tstats,10000)##modules with vertex in EntrezID (FIS), graph to be test and meth stats (P6),number of randomizations

stats.perm<-as.data.frame(Stats$fdr)
##now get correct direction

##for each component of the module, extract meth val and locate in its module

stats.sp<-list()
##get the genes in P6 set
stats.sp<-lapply(lv, function(x) {
  
  return(intEpi.o$statM[x,"t"])
  return(rownames(intEpi.o$statM[x]))} )

#re-scale



library(scales)
stats.sp2<-lapply(stats.sp, function(x) {return(rescale(x, to = c(-2, 2)))} )



#get original stats 
ori.stats<-lapply(Spain.Epi$topmod, function(x)return(x$'stat(DNAm)'))

##re-scale
ori.stats2<-lapply(ori.stats, function(x) return(rescale(x,to=c(-2,2))))



##Toni's option

##Sobre replicació moduls: de cada modul fer un vector dels estadístics dels gens (crec que era un tvalue: te significació i sentit)
##i comparar amb un paired t-test el vector del dataset 1 vs el 2. 
#No hauria de donar significatiu per confirmar que estem veient el mateix tipus de mòdul.

##wilcoxon rank test

#myfunction <- function(x,y) {
 # test = wilcox.test(x, y, paired=T)
  #return(test$p.value)

#}

## t test
myfunction <- function(x,y) {
  test = t.test(x, y, paired=T)
  return(test$p.value)
}


results.wilcox<-mapply(myfunction, x = ori.stats2, y = stats.sp2, SIMPLIFY = FALSE)

results.wilcox<-as.data.frame(t(as.data.frame(results.wilcox)))


Stats.r<-merge(stats.perm, results.wilcox, by="row.names")
colnames(Stats.r)<-c("Module","pval_Permut","pval_Direction")


write.csv(Stats.r,"/mnt/md127/Analysis_Projects/FIS_Meth/Results_antiTNF/delta/Validation_Modules_delta_P6.csv",row.names = F)
#write.csv(Stats.r,"/mnt/md127/Analysis_Projects/FIS_Meth/Results_pre_selected/Quant_baseline_all_treatments/validated_modules/Tats_validation_Quant_Baseline_all_treatments.csv",row.names = F)


#######################Perm graph

pdf("Perm_Plot_validation_DeltaDAS28_antiTNF.pdf")

for( i in 1:length(Stats$Random)){

xlims.0<-min(Stats$Random[[i]])
xlims.1<-Stats$Observed[[i]]+5

h<-hist(Stats$Random[[i]], col="lightblue",
     xlab = "FDR estimate",
     xlim = c(xlims.0,xlims.1),
     main = paste("Validation DeltaDAS28 wk0 antiTNF ",names(Stats$fdr)[i], sep="")
)

abline(v = Stats$Observed[[i]], col = "red", lwd = 4)
abline(v = mean(Stats$Random[[i]]), col = "black", lwd = 4)


xfit<-seq(min(Stats$Random[[i]]),max(Stats$Random[[i]]),length=40)
yfit<-dnorm(xfit,mean=mean(Stats$Random[[i]]),sd=sd(Stats$Random[[i]]))
yfit <- yfit*diff(h$mids[1:2])*length(Stats$Random[[i]])

lines(xfit, yfit, col="black", type="l",lwd=2,lty=2) 
arrows(mean(Stats$Random[[i]]),max(h$counts)-mean(h$counts),Stats$Observed[[i]],max(h$counts)-mean(h$counts),code=1)
arrows(mean(Stats$Random[[i]]),max(h$counts)-mean(h$counts),Stats$Observed[[i]],max(h$counts)-mean(h$counts),code=2)


quantiles <- quantile(Stats$Random[[i]], prob=0.95)
abline(v = as.numeric(quantiles), col = "green", lwd = 4)


}
dev.off()

####

#boxplot(c(ori.stats2[[2]],stats.sp2[[2]])~c(rep("A",length(ori.stats2[[2]])),rep("B",length(ori.stats2[[2]]))))


##try to plot Modules GOs in Validation set graph


##plots
library(FEM)
library(marray)
library(corrplot)

setwd("/mnt/md127/Analysis_Projects/FIS_Meth/Results_antiTNF/delta/Vals/")##prior to plot!!!


selected.mods<-Stats.r[Stats.r$pval_Permut<.05 & Stats.r$pval_Direction>0.05,"Module"]
sels<-which(names(Spain.Epi$topmod) %in% selected.mods)


##plot selected modules in validation set



#Spain.Epi$topmod[[1]]

#EntrezID   Symbol  stat(DNAm)      P(DNAm)  stat(Int)
#2475       2475     MTOR -3.83139428 0.0002889959 3.83139428
#7248       7248     TSC1  0.53949448 0.5913830722 0.53949448
#1978       1978 EIF4EBP1 -0.16312998 0.8709205750 0.16312998

###build a list with the modules selected

Modules.to.test<-data.frame(EntrezID=lv[[1]],"stat(DNAm)"=as.numeric(stats.sp[[1]]))

Modules.to.test<-Map(cbind, "EntrezID"=lv, "stat(DNAm)" = stats.sp)
Modules.to.test<-lapply(Modules.to.test, function(x) {return(as.data.frame(x))} )


##merge the two lists of data frames by EntrezID
library(purrr)
library(dplyr)

List<-map2(Modules.to.test, Spain.Epi$topmod, inner_join, by = "EntrezID")

List<-lapply(List, function(x) setNames(x, sub("\\.x", "", names(x))))
List<-lapply(List, function(x){
  
            rownames(x)<-x$EntrezID
            return(x)
  } )


##select the ones that want to plot
for(m in 1:length(sels)){
  FemModShow.val(List[[sels[m]]],
             name=names(List)[sels[m]],EpiMod.o,mode="Epi")
  
}



#Modules.to.test<-data.frame(EntrezID=lv[[3]],"stat(DNAm)"=as.numeric(stats.sp[[3]]))
#colnames(Modules.to.test)[2]<-"stat(DNAm)"
#Modules.to.test$EntrezID<-as.factor(Modules.to.test$EntrezID)

#Modules.to.test<-merge(Modules.to.test,Spain.Epi$topmod[[3]][,c("EntrezID","Symbol","P(DNAm)","stat(Int)")],by="EntrezID")
#Modules.to.test$EntrezID<-as.factor(Modules.to.test$EntrezID)
#rownames(Modules.to.test)<-Modules.to.test$EntrezID

#source("/mnt/md127/Analysis_Projects/FIS_Meth/Scripts_last/FemModShow_val.R")
#FemModShow.val(Modules.to.test,name="Test2",EpiMod.o,mode="Epi")






##GO enrichment for each module

#####GO enrichment analysis
library(clusterProfiler)
library(org.Hs.eg.db)
#library(GOstats)
library(GO.db)
library(DOSE)
library(annotate)


#plot.list<-list()


for(m in 1:length(sels)){
  
  
  genesid<-as.character(Spain.Epi$topmod[[sels[m]]][,"Symbol"])
  
  genesid<-genesid[ genesid != "" ]##remove empty elements from a vector
  genesid <- genesid[!is.na(genesid)]
  
  eg<-bitr(genesid, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  
  ego2 <- enrichGO(gene         = eg$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENTREZID',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable=T
  )
  
  png(paste0("GO_Module","_","Delta_w0_Vals","_",names(Spain.Epi$topmod)[sels[m]],".png",sep=""),width = 1200, height = 800)
  
  p<-dotplot(ego2, title=names(Spain.Epi$topmod)[sels[[m]]],showCategory=35)
  
  print(p)
  
  dev.off()
  # myplots[[i]] <- p
}




