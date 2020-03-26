library(Seurat)
# take normalised data
count_norm <- GetAssayData(seurat_object)
write.table(count_norm, 'cellphonedb_count.txt', sep='\t', quote=F)
# generating meta file
meta_data <- cbind(rownames(seurat_object@meta.data), cluster=seurat_object@meta.data[,'Cluster', drop=F])   #####  cluster is the user’s specific cluster column
write.table(meta_data, 'cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)

means<-read.table("./cellphonedb-result/means.txt",sep = "\t",header=T,check.names = F)
means_sim<-means[,12:60]


pvalues<-read.table("./cellphonedb-result/pvalues.txt",sep = "\t",header=T,check.names = F)
pvalues<-pvalues[match(means$id_cp_interaction,pvalues$id_cp_interaction),]
pvalues_sim<-pvalues[,12:60]


cluster<-unique(unlist(strsplit(colnames(means_sim),split = "\\|"))) 
l<-length(cluster)

means_sel<-matrix(nrow = nrow(means_sim),ncol = ncol(means_sim))
for(i in 1:nrow(means_sim)){
  for(j in 1:ncol(means_sim)){
    ifelse(pvalues_sim[i,j]<=0.05,means_sel[i,j]<-means_sim[i,j],means_sel[i,j]<-0)
  }
}
colnames(means_sel)<-colnames(means_sim)

row_sum<-rowSums(means_sel)
index<-head(order(row_sum,decreasing = T),5)
means$interacting_pair[index]
means_i<-means[index,]
selected_rows<-as.character(means_i$interacting_pair)

pvalues_sim_logi<-pvalues_sim<=0.05
pvalues_row_sum<-rowSums(pvalues_sim_logi)
index<-names(head(sort(pvalues_row_sum,decreasing = T),5))
pvalues_i<-pvalues[index,]
selected_rows<-pvalues_i$interacting_pair
selected_rows<-as.character(selected_rows)