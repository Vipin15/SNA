library(igraph)

L1_comm=read.csv(file='/disk1/Workspace/SNA_maxis_oct/KL_Comm_details.csv',header=T)
CD_L1_comm=subset(L1_comm,!grepl("_",L1_comm$Comm_label))

Inf_CD_L1=data.frame()


for(i in 1:length(CD_L1_comm$Comm_label)){
  
  C1=induced.subgraph(G2, which(membership(Fast_community_unfolding)==as.numeric(as.character(CD_L1_comm$Comm_label[i]))))
  
  V(C1)$Comm_label=as.numeric(as.character(CD_L1_comm$Comm_label[i]))
  
  
  V(C1)$Degree=degree(C1, v=V(C1),mode=c("all"),loops=F)
  
  #CALCULATING METRICS
  
  #Calculating weighted betweeness vector centrality for each node
  V(C1)$betweennes = estimate_betweenness(C1,vids = V(C1),directed = F,weights = E(C1)$weight,cutoff = 6)
  V(C1)$betweennes=replace(V(C1)$betweennes,is.nan(V(C1)$betweennes),0)
  V(C1)$betweennes=round(1+(V(C1)$betweennes - min(V(C1)$betweennes))*9/diff(range(V(C1)$betweennes)),3)
  V(C1)$betweennes=replace(V(C1)$betweennes,is.nan(V(C1)$betweennes),0)
  
  #Calculating weighted eigen_centrality vector centrality for each node
  EV=eigen_centrality(C1,weights = E(C1)$weight,scale = F)
  EV=as.data.frame(EV$vector)
  V(C1)$EV=EV$`EV$vector`
  rm(EV)
  V(C1)$EV=replace(V(C1)$EV,is.nan(V(C1)$EV),0)
  V(C1)$EV=round(1+(V(C1)$EV - min(V(C1)$EV))*9/diff(range(V(C1)$EV)),3)
  V(C1)$EV=replace(V(C1)$EV,is.nan(V(C1)$EV),0)
  
  
  #Clustering coefficent of each node
  V(C1)$Clustering_coef=transitivity(C1,type = "local")
  V(C1)$Clustering_coef=replace(V(C1)$Clustering_coef,is.nan(V(C1)$Clustering_coef),0)
  V(C1)$Clustering_coef=round(1+(V(C1)$Clustering_coef - min(V(C1)$Clustering_coef))*9/diff(range(V(C1)$Clustering_coef)),3)
  V(C1)$Clustering_coef=replace(V(C1)$Clustering_coef,is.nan(V(C1)$Clustering_coef),0)
  
  
  #K_core
  V(C1)$K_core=coreness(C1,mode = "all")
  V(C1)$K_core=replace(V(C1)$K_core,is.nan(V(C1)$K_core),0)
  V(C1)$K_core=round(1+(V(C1)$K_core - min(V(C1)$K_core))*9/diff(range(V(C1)$K_core)),3)
  V(C1)$K_core=replace(V(C1)$K_core,is.nan(V(C1)$K_core),0)
  
  #Wt. Degree
  V(C1)$Wt.Degree=graph.strength(C1,vids = V(C1),mode = c("all"),loops = F)
  V(C1)$Wt.Degree=replace(V(C1)$Wt.Degree,is.nan(V(C1)$Wt.Degree),0)
  V(C1)$Wt.Degree=round(1+(V(C1)$Wt.Degree - min(V(C1)$Wt.Degree))*9/diff(range(V(C1)$Wt.Degree)),3)
  V(C1)$Wt.Degree=replace(V(C1)$Wt.Degree,is.nan(V(C1)$Wt.Degree),0)
  
  
  #Wt. Cloeseness
  V(C1)$Closenes=closeness(C1,vids = V(C1),mode = c("all"),weights = E(C1)$weight)
  V(C1)$Closenes=replace(V(C1)$Closenes,is.nan(V(C1)$Closenes),0)
  V(C1)$Closenes=round(1+(V(C1)$Closenes - min(V(C1)$Closenes))*9/diff(range(V(C1)$Closenes)),3)
  V(C1)$Closenes=replace(V(C1)$Closenes,is.nan(V(C1)$Closenes),0)
  
  V(C1)$Score=V(C1)$K_core+V(C1)$Clustering_coef+V(C1)$Wt.Degree+V(C1)$Closenes+V(C1)$EV+V(C1)$betweennes
  V(C1)$Score=round(1+(V(C1)$Score - min(V(C1)$Score))*99/diff(range(V(C1)$Score)),3)
  V(C1)$Score=replace(V(C1)$Score,is.nan(V(C1)$Score),0)
  
  
  library(dplyr)
  detach(package:dplyr)
  
  #Nodelist
  nodes_influencers=get.data.frame(C1, what="vertices")
  library(dplyr)
  nodes_influencers = nodes_influencers %>% select(name,
                                                   Degree,
                                                   betweennes,
                                                   EV,
                                                   Clustering_coef,
                                                   K_core,
                                                   Wt.Degree,
                                                   Closenes,
                                                   Score)
  
  
  Inf_CD_L1 = rbind(Inf_CD_L1, nodes_influencers)
  
  
  print("CD_L1_percentage")
  print(nrow(Inf_CD_L1)*100/sum(CD_L1_comm$Size))
  
  
}
Inf_CD_L1=aggregate(.~name,data=Inf_CD_L1,median)
