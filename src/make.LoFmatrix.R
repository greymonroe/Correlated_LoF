#### make LoFmarix ###

load("RData/gene_model_stops")

cl=makeCluster(7)
registerDoSNOW(cl)
clusterExport(cl, list("gene_model_stops"))

gene_model_stops_ecotype_count<-foreach(g=1:length(gene_model_stops)) %dopar% {
  
  return(length(gene_model_stops[[g]]))
  
}

stopCluster(cl)

table(unlist(gene_model_stops_ecotype_count))
which(unlist(gene_model_stops_ecotype_count)==331)


gene_model_stops_df<-do.call("cbind", gene_model_stops[-which(unlist(gene_model_stops_ecotype_count)==331)])
gene_model_stops_df<-data.frame(gene_model_stops_df)
colnames(gene_model_stops_df)<-genes[-which(unlist(gene_model_stops_ecotype_count)==331)]

gene_model_maxes<-apply(gene_model_stops_df, 2, function(x){
  return(x<max(x, na.rm=T)*.9)
})

gene_model_maxes<-data.frame(gene_model_maxes)

genes_plain<-unique(gsub("(*)\\..+", "\\1", genes))


cl=makeCluster(7)
registerDoSNOW(cl)
clusterExport(cl, list("gene_model_maxes","genes_plain"))

genes_plain_lof<-foreach(g=1:length(genes_plain), .combine=cbind) %dopar% {
  
  conc<-apply(gene_model_maxes[,grep(genes_plain[g], colnames(gene_model_maxes)), drop=F], 1, function(x) T %in% x)
  return(conc)
  
}

stopCluster(cl)

LoFmatrix<-data.frame(genes_plain_lof)
colnames(LoFmatrix)<-genes_plain

save(LoFmatrix, file="RData/LoFmatrix")
