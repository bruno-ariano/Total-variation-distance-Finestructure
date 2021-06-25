###############################################################################
#Bruno Ariano                                                                 #
#Total variation distance(TVD) using Chromopainter and Finestructure results  #
#mail: bruno.ariano.87@gmail.com; arianob@tcd.ie                              #
###############################################################################

library(dplyr)
library(readr)
library(ape)
library(phytools)
library(geiger)
#This library was taken from the Finestructure website
source("/home/bruno/Genomics_folder/Software/fs_4.1.1/FinestructureLibrary.R")

finestr = "/home/bruno/Genomics_folder/Chromopainter/last_run/n563_M9.50e05.chunkcount.x100k_y100k.mcmc.T.k2.tree"
treexml<-xmlTreeParse(finestr)
ttree<-extractTree(treexml)

chunklength = "/home/bruno/Genomics_folder/Chromopainter/last_run/n563_M9.50e05_combine.chunklengths.out"
data_chunks = read.table(chunklength, header = T, row.names = 1)

#Normalization step
data_chunks_norm = t(apply(data_chunks, 1, function(x) x/sum(x)))


###### To define groups I cut the tree at arbitrary poit using the treeSlice function
tree_sliced = treeSlice(rescale(ttree, "delta", 5),prompt = T,orientation="rootwards")
root_nodes = as.numeric(tree_sliced$tip.label)

#In case the above code give error use the below (slice indicate where you want to cut the tree)
#tree_sliced = treeSlice(ttree,prompt = F,slice = 120000) 
# root_nodes = c()
# for(record in tree_sliced){
#   root_nodes = c(root_nodes, findMRCA(ttree, tip = record$tip.label))
# }
                          
#root_nodes = c()
#for(record in tree_sliced){
#  root_nodes = c(root_nodes, findMRCA(ttree, tip = record$tip.label))
#}

#Here I store all pairs of combination between nodes(groups) of the tree
combn_nodes = combn(root_nodes, 2)

#Function to calculate the copying vectors
make_copy_vector = function(matrix_chunk, label_group1, label_group2, root_nodes){
  #First I extract the matrix with the labels for each group
  group_matrix = matrix_chunk[label_group1,]
  #First I average by row
  group_matrix_mean = apply(group_matrix,2,mean)
  #Then I sum all columns(Donors)
  copy_vector1 = c()
  for(node in root_nodes){
      tree = extract.clade(ttree,node)
      copy_vector1 = c(copy_vector1, sum(group_matrix_mean[tree$tip.label]))
  }
  group_matrix2 = matrix_chunk[label_group2,]
  group_matrix_mean2 = apply(group_matrix2,2,mean)
  copy_vector2 = c()
  for(node in root_nodes){
    tree = extract.clade(ttree,node)
    copy_vector2 = c(copy_vector2, sum(group_matrix_mean2[tree$tip.label]))
  }
  return(list(copy_vector1,copy_vector2))
}
                           
#All TVD values will be stored in a square matrix with row and col. names given by node label
tvd_matrix = matrix(rep(0,length(root_nodes)**2),nrow = length(root_nodes))
colnames(tvd_matrix) = root_nodes
rownames(tvd_matrix) = root_nodes
                        
                           
for(i in seq(ncol(combn_nodes))){
  
  tree1 = extract.clade(ttree,combn_nodes[1,i])
  tree2 = extract.clade(ttree,combn_nodes[2,i])
  copy_vect = make_copy_vector(data_chunks_norm, tree1$tip.label, tree2$tip.label)
  tvd = 0.5 * sum(abs(copy_vect[[1]] - copy_vect[[2]]))
  
  #I create 200 random groups by sampling from the two selected groups. Depending on the clusters size you might want to increase the permutation number.
  tvd_random = c()
  all_label = c(tree1$tip.label, tree2$tip.label)
  for(j in seq(200)){
    group1_label_rand = sample(all_label, length(tree1$tip.label))
    group2_label_rand = all_label[!all_label %in% group1_label_rand]
    copy_vect_rand = make_copy_vector(data_chunks_norm, group1_label_rand, group2_label_rand)
    tvd_random = c(tvd_random, 0.5 * sum(abs(copy_vect_rand[[1]] - copy_vect_rand[[2]])))
  }
  wt = wilcox.test(tvd_random,mu = tvd)
  print(wt)
  tvd_matrix[as.character(combn_nodes[1,i]),as.character(combn_nodes[2,i])] =tvd 
  tvd_matrix[as.character(combn_nodes[2,i]),as.character(combn_nodes[1,i])] =tvd 
}

                           
