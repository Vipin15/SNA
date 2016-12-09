library(igraph)

L1_comm=read.csv(file='/disk1/Workspace/SNA_maxis_oct/KL_Comm_details.csv',header=T)

CD_L2_comm=subset(L1_comm,grepl("_",L1_comm$Comm_label))

df=as.data.frame(unique(sub("_.*", "", CD_L2_comm$Comm_label)))
names(df)="id"
Influencers_subcomm = data.frame()

for (i in df$id) {
     C2=induced.subgraph(G2, which(membership(Fast_community_unfolding)==i))
       CD_2=cluster_louvain(C2,weights = E(C2)$weight)

         q = data.frame(
             factor(),
                 integer(),
                     integer(),
                         integer(),
                             integer(),
                                 integer(),
                                     integer(),
                                         integer(),
                                             factor()
                                               )
         df=as.data.frame(sizes(CD_2))

         df=subset(df,df$Freq>2)

         names(df)=c("Comm_label","Size")


for (j in df$Comm_label) {

     C1=induced.subgraph(C2, which(membership(CD_2)==j))

      V(C1)$Comm_label=paste0(i,sep="_",j)


          V(C1)$Degree=degree(C1, v=V(C1),mode=c("all"),loops=F)

              #CALCULATING METRICS

          
#Calculating weighted betweeness vector centrality for each node
              V(C1)$betweennes = estimate_betweenness(C1,vids = V(C1),directed = F,weights = E(C1)$weight,cutoff = 4)
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


 V(C1)$Score=V(C1)$K_core+V(C1)$Clustering_coef+V(C1)$Wt.Degree+V(C1)$EV+V(C1)$betweennes+V(C1)$Closenes
     V(C1)$Score=round(1+(V(C1)$Score - min(V(C1)$Score))*99/diff(range(V(C1)$Score)),3)
         V(C1)$Score=replace(V(C1)$Score,is.nan(V(C1)$Score),0)



 #Nodelist
library(dplyr)
  detach(package:dplyr)


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

  
    q = rbind(q, nodes_influencers)
      
        print("SubComm_perc")
          print(length(unique(q$Comm_label))*100/length(unique(membership(CD_2))))
            
              
}

Influencers_subcomm = rbind(Influencers_subcomm, q)
rm(q)

print("perc_nodes_covered")
print(nrow(Influencers_subcomm)*100/sum(CD_L2_comm$Size))
}
