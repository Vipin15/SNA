#Influencers subcommunity
library(igraph)
library(dplyr)
start.time <- Sys.time()

#closely_knit_labels=subset(Non_isolated_comm,grepl("Closely",Non_isolated_comm$Status)|grepl("Small",Non_isolated_comm$Status)|grepl("Medium",Non_isolated_comm$Status))

#a=subset(isolated_comm,isolated_comm$Size>2)

#Cycle1=subset(Non_isolated_comm,Non_isolated_comm$Status=="Large & Loosely Knit"|Non_isolated_comm$Status=="VERY LOW DENSITY")

L1_comm=read.csv(file='/disk1/Workspace/SNA_maxis_oct/Sabah_Comm_details.csv',header=T)


closely_knit_labels=subset(L1_comm,L1_comm$Secluded=='NO'&!grepl("_",L1_comm$Comm_label))

a=subset(L1_comm,L1_comm$Secluded=='YES')

Cycle1=subset(L1_comm,grepl("_",L1_comm$Comm_label))

Cycle1=Cycle1[order(-Cycle1$Size,-Cycle1$Density),]

df=as.data.frame(unique(sub("_.*", "", Cycle1$Comm_label)))
names(df)="id"

total_nodes=sum(closely_knit_labels$Size,a$Size,Cycle1$Size)

Total_influencers=25000

influencers_subc=ceiling((Total_influencers*sum(Cycle1$Size))/total_nodes)


Influencers_subcomm = data.frame(
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

for (i in df$id) {
  
  if (sum(Influencers_subcomm$Influencer == "YES")>=influencers_subc){
    next()
  }
  if (sum(Influencers_subcomm$Influencer=="YES"|Influencers_subcomm$Influencer=="DONE")/sum(Cycle1$Size)>0.7){
    next()
  }
  
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

df=df[order(-df$Size),]

  for (j in df$Comm_label) {
    
    if (sum(q$Influencer == "YES")>=influencers_subc){
      next()
    }
    C1=induced.subgraph(C2, which(membership(CD_2)==j))
    
    #if (sum(Influencers_subcomm$Influencer=="YES"|Influencers_subcomm$Influencer=="DONE")/sum(Cycle1$Size)>0.7){
    #next()
    # }
   # V(C1)$Original.Comm_label=i
    #V(C1)$SubComm_label=j

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
    
    
    #Edgelist
    
    detach(package:dplyr)
    Edges_inf=as_data_frame(C1, what="edges")
    edge_weights=quantile(Edges_inf$weight)
    
    #Nodelist
    
    nodes_influencers=as_data_frame(C1, what="vertices")
    
    
library(dplyr)
 nodes_influencers=nodes_influencers %>% select(name,Degree,Comm_label,Score)


     nodes_influencers$Influencer="NO"
    
    
    
    

    community <- nodes_influencers %>% 
      arrange(desc(Score)) %>% 
      filter(Degree > 1) %>% 
      select(name)
    
    
    
    threshold=ceiling(influencers_subc*nrow(nodes_influencers)/sum(Cycle1$Size))
    
    for (person in community$name) {
      
      if ( nodes_influencers[nodes_influencers$name == person, 'Influencer'] == 'DONE'|
           nodes_influencers[nodes_influencers$name == person, 'Influencer'] == 'REJECTED'|
           sum(nodes_influencers$Influencer=="YES")>=threshold){
        next()
      }
      
      nodes_influencers[nodes_influencers$name == person, 'Influencer'] <- 'YES'
      
      subdf <- Edges_inf %>%
        filter(from == person | to == person) %>%
        filter(weight >= edge_weights[3])
      
      # neighbors <- union(subdf$from, subdf$to)
      
      neighbors <- setdiff(union(subdf$from, subdf$to), person)
      
      if (length(neighbors) == 0) {
        
        nodes_influencers[nodes_influencers$name == person, 'Influencer'] <- 'REJECTED'
        
      } else {
        
        overlap <- (length(neighbors) - 
                      sum(nodes_influencers[nodes_influencers$name %in% neighbors, 'Influencer'] == 'NO' |
                            nodes_influencers[nodes_influencers$name %in% neighbors, 'Influencer'] == 'REJECTED')) / length(neighbors)
        
        if (overlap <= 0.5) {
          nodes_influencers[nodes_influencers$name == person, 'Influencer'] <- 'YES'
        } else {
          nodes_influencers[nodes_influencers$name == person, 'Influencer'] <- 'REJECTED'
        }
        
        nodes_influencers[(nodes_influencers$name %in% neighbors & 
                             (nodes_influencers$Influencer == 'NO' |
                                nodes_influencers$Influencer == 'REJECTED')), 'Influencer'] <- 'DONE'
      }
      
    }
   
 nodes_influencers=nodes_influencers %>% select(name,Degree,Comm_label,Score,Influencer)

    q = rbind(q, nodes_influencers)
    
   print("SubComm_percentage")
   print(length(unique(q$Comm_label))*100/length(unique(membership(CD_2))))
    
  }
  
  Influencers_subcomm = rbind(Influencers_subcomm, q)
  rm(q)
  print("Subcomm_Comm_percentage")
  print(length(unique(Influencers_subcomm$Comm_label))*100/nrow(df))
  
  print("Subcomm_perc_influencers_covered")
  print(sum(Influencers_subcomm$Influencer=="YES")*100/influencers_subc)
}


#rm(C1,C2,CD_2,i,perc_influencers_covered,subcomm_perc,j,overlap,nodes_influencers,neighbors,subdf,person,Edges_inf,edge_weights,community)

end.time <- Sys.time()

subcomm.time=round(difftime(end.time,start.time),2)
cat("subcomm.time=",subcomm.time,"\n")

rm(end.time,start.time)
