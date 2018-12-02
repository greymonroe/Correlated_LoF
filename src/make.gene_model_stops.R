### make gene_model_stops ###

path_to_CDS_fastas="~/Dropbox/Research/LOF_questions/Data/final_CDS/"
CDS_files<-list.files(path_to_CDS_fastas, full.names = T)
genes<-gsub("_CDS.fasta", "", list.files(path_to_CDS_fastas, full.names = F))

path_to_p3_fastas="~/Dropbox/Research/LOF_questions/Data/final_3pUTR/"
p3_files<-list.files(path_to_p3_fastas, full.names = T)

# scan fasta files to identify gene_model_stops --------------------------------------------------------

cl=makeCluster(7)
registerDoSNOW(cl)
clusterExport(cl, list("CDS_files","p3_files", "genes"))
clusterEvalQ(cl, library(seqinr))

gene_model_stops<-foreach(g=1:length(genes)) %dopar% {
  
  gene<-genes[g]
  
  print(g)
  CDS_file=paste0(path_to_CDS_fastas, gene, "_CDS.fasta")
  p3_file=paste0(path_to_p3_fastas, gene, "_3pUTR.fasta")
  
  CDS_fasta<-read.fasta(CDS_file)
  
  if (p3_file %in% p3_files){
    p3_fasta<-read.fasta(p3_file)
    p3<-T
  } else p3=F
  
  
  ends<-c()
  for(e in 1:length(CDS_fasta)){
    ecotype<-names(CDS_fasta)[e]
    prot<-translate(CDS_fasta[[e]])
    
    if (p3){
      p3_dna<-p3_fasta[[e]]
      if (length(p3_dna)>3){
        prot<-c(prot, translate(p3_dna))
      }
    }
    
    end=which(prot=="*")[1]
    #cat( "\t ecotype #", e, "...", ecotype, "... gene #", g , "...", gene_model, "\n")
    
    ends<-c(ends, end)
  }
  
  return(ends)
  
}

stopCluster(cl)

#uncomment to save new gene_model_stops
save(gene_model_stops, file="RData/gene_model_stops")
