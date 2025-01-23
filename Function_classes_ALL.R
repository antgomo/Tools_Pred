
###FUNCTION DEFINITION

####Select classes in ALL set, avoid those classes with less than 6 samples selecting the previous K in M3C
####Using K in Set1, Mirror in Test set, Set2, avoid classes with less than 6 samples, removing those samples in these low representative
####classes and running M3C without these outliers if we are delaing with no more than 4 outliers, othewise, we mark this combinatorial analysis 
####approach as not valid

####

#hi.adj discovery set adjusted data, aka Set1
#lo.adj test set adjusted data, aka Set2
#A.set1 discovery set samples
#A.set2 test set samples

#imid imid tested
#genes.sel.hi genes used to test in discovery : DEG, VAR,... according combinations. Keep in mind that when DEG is used,
##we are using those DEG found in discovery ser independently from DEG in test
#genes.sel.lo genes used to test in test set : DEG, VAR,... according combinations. Keep in mind that when DEG is used,
##we are using those DEG found in discovery independently from DEG in test


class.pipeline<-function(hi.adj,A,imid,genes.sel.hi){
  
  
  ###Discovery set data
  
  data.tp<-hi.adj[rownames(hi.adj) %in% genes.sel.hi,colnames(hi.adj) %in% A]
  
  ##run M3C
  res.m3c <- M3C(data.tp, cores=7,seed = 123, removeplots = F)
  stats<-res.m3c$scores #store stats
  
  selected.stat<-stats[stats$RCSI==max(stats$RCSI),c("K","RCSI","MONTECARLO_P")]
  pval<-round(selected.stat$MONTECARLO_P, digits = 3)
  
  j<-selected.stat$K
  sel.c<-selected.stat$K## This will be K selected and uset in Set2 to fix, aka Mirroring

    data <- res.m3c$realdataresults[[j]]$ordered_data # this is the data
  annon <- res.m3c$realdataresults[[j]]$ordered_annotation # this is the annotation
  #write annot
  write.csv(annon,paste0(output_dir,"/","Annon_classes_ALL_original",my.imids[r],"_",my.options[i],".csv",sep=""))  
  write.csv(stats,paste0(output_dir,"/","Stats_classes_ALL",imid,"_",my.options[i],".csv",sep=""))
  
  #####avoid problems with low classes number
  
  ##avoid classes less than 6 samples
  if(min(table(annon))<6){
    
    ##while less than 6 samples in one class, go down to select another K to fix
    #end when more than 5 samples in one class
    while(min(table(annon))<6){

      ##go for K-1
      selected.stat<-stats[stats$K==j-1,c("K","RCSI","MONTECARLO_P")]
      pval<-round(selected.stat$MONTECARLO_P, digits = 3)
      
      j<-selected.stat$K
      ##new sel.c to mirror in Set2
      sel.c<-selected.stat$K

      data <- res.m3c$realdataresults[[j]]$ordered_data # this is the data
      annon <- res.m3c$realdataresults[[j]]$ordered_annotation # this is the annotation
      ccmatrix <- res.m3c$realdataresults[[j]]$consensus_matrix # this is the consensus matrix
      # normalise and scale the data
      data <- t(scale(t(data))) # z-score normalise each row (feature)
      data <- apply(data, 2, function(x) ifelse(x > 4, 4, x)) # compress data within range
      data <- apply(data, 2, function(x) ifelse(x < -4, -4, x)) # compress data within range
      
      write.csv(annon,paste0(output_dir,"/","Annon_classes_ALL",my.imids[r],"_",my.options[i],".csv",sep=""))  
      write.csv(stats,paste0(output_dir,"/","Stats_classes_ALL",imid,"_",my.options[i],".csv",sep=""))
      
      cons<-ggsci::pal_futurama()(max(levels(annon$consensuscluster)))[as.factor(annon$consensuscluster)]#5 clusters
      
      
      spcol <- cons
      
      cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
      png(paste0(output_dir,"/","heatmap_Classes_ALL_samples",imid,"_",my.options[i],".png",sep=""),width = 800, height = 800)
      
      
      heatmap.2x(data, col=rev(cols), Colv = NA,Rowv=TRUE, scale="none", ColSideColors=spcol,
                 trace="none", dendrogram="both", 
                 cexRow=1, cexCol=1.4,
                 main=paste(imid,"ALL samples genes",my.options[i],pval,sep=" "),
                 labCol=NA,labRow=NA, density.info="none",
                 hclust=function(x) hclust(x,method="complete"),
                 distfun=function(x) as.dist((1-cor(t(x)))/2)
                 
      )
      dev.off()
    }
  }else{
    ###If orginal K in M3C, has more than 5 samples in one class, avoid first loop
    # normalise and scale the data
    data <- t(scale(t(data))) # z-score normalise each row (feature)
    data <- apply(data, 2, function(x) ifelse(x > 4, 4, x)) # compress data within range
    data <- apply(data, 2, function(x) ifelse(x < -4, -4, x)) # compress data within range
    
    write.csv(stats,paste0(output_dir,"/","Stats_classes_ALL",imid,"_",my.options[i],".csv",sep=""))
    write.csv(annon,paste0(output_dir,"/","Annon_classes_ALL",my.imids[r],"_",my.options[i],".csv",sep=""))  
    
    cons<-ggsci::pal_futurama()(max(levels(annon$consensuscluster)))[as.factor(annon$consensuscluster)]#5 clusters
    
    
    spcol <- cons
    
    cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
    png(paste0(output_dir,"/","heatmap_classes_ALL_samples",imid,"_",my.options[i],".png",sep=""),width = 800, height = 800)
    
    
    heatmap.2x(data, col=rev(cols), Colv = NA,Rowv=TRUE, scale="none", ColSideColors=spcol,
               trace="none", dendrogram="both", 
               cexRow=1, cexCol=1.4,
               main=paste(imid,"ALL samples",my.options[i],pval,sep=" "),
               labCol=NA,labRow=NA, density.info="none",
               hclust=function(x) hclust(x,method="complete"),
               distfun=function(x) as.dist((1-cor(t(x)))/2)
               
    )
    dev.off()
  }
  rm(res.m3c)
  rm(stats)
  rm(selected.stat)
  rm(data)
  gc() 
  
  
  
###end of process
}
