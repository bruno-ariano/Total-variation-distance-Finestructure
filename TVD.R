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
source("/home/bruno/Genomics_folder/Software/fs_4.1.1/FinestructureLibrary.R")

finestr = "/home/bruno/Genomics_folder/Chromopainter/last_run/n563_M9.50e05.chunkcount.x100k_y100k.mcmc.T.k2.tree"
treexml<-xmlTreeParse(finestr) ## read the tree as xml format
ttree<-extractTree(treexml) ## extract the tree into ape's phylo format

chunklength = "/home/bruno/Genomics_folder/Chromopainter/last_run/n563_M9.50e05_combine.chunkcounts.out"
data_chunks = read.table(chunklength,header = T, row.names = 1)

#Normalization step
data_chunks_norm = t(apply(data_chunks, 1, function(x) x/sum(x)))


###### To define groups I cut the tree at arbitrary poit using the treeSlice function
tree_sliced = treeSlice(ttree,prompt = T)

root_nodes = c()
for(record in tree_sliced){
  root_nodes = c(root_nodes, findMRCA(ttree, tip = record$tip.label))
}

#Here I store all pairs of combination between nodes(groups) of the tree
combn_nodes = combn(root_nodes, 2)

#two vectors will store the real groups and the random generated ones
list_group = c()
list_group_random = c()

for(i in seq(ncol(combn_nodes))){
  tree1 = extract.clade(ttree,combn_nodes[1,i])
  group_matrix1 = data_chunks_norm[tree1$tip.label,]
  copy_vect1 = apply(group_matrix1,2,mean)
  
  tree2 = extract.clade(ttree,combn_nodes[2,i])
  group_matrix2 = data_chunks_norm[tree2$tip.label,]
  copy_vect2 = apply(group_matrix2,2,mean)
  
  tvd = 0.5 * sum(abs(copy_vect1 - copy_vect2))
  
  #I create 200 random groups by sampling from the two above ones
  tvd_random = c()
  all_label = c(tree1$tip.label, tree2$tip.label)
  for(j in seq(200)){
    group1_label_rand = sample(all_label, length(tree1$tip.label))
    group2_label_rand = all_label[!all_label %in% group1_label_rand]
    group_rand_matrix1 = data_chunks_norm[group1_label_rand,]
    group_rand_matrix2 = data_chunks_norm[group2_label_rand,]
    copy_rand_vect1 = apply(group_rand_matrix1,2,mean)
    copy_rand_vect2 = apply(group_rand_matrix2,2,mean)
    tvd_random = c(tvd_random, 0.5 * sum(abs(copy_rand_vect1 - copy_rand_vect2)))
    
  }
  wt = wilcox.test(tvd_random,mu = tvd)
  print(wt)
}
