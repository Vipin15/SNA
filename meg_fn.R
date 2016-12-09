
library(igraph)
library(dplyr)
options(scipen = 999)


Node_Comm_fn <- function(G,Comm_obj)
{
  
  
  df=as.data.frame(sizes(Comm_obj))
  
  
  names(df)=c("Comm_label","Size")
  
  
  a=data.frame(
    numeric()
  )
  
  for(j in 1:nrow(df)){
    C=induced.subgraph(G, which(membership(Comm_obj)==j))
    mean_deg_bfr=mean(V(C)$Degree)
    V(C)$Degree=degree(C, v=V(C),mode=c("all"),loops=F)
    Secluded=ifelse((mean_deg_bfr-mean(V(C)$Degree))==0,"YES","NO")
    
    p = graph.density(C,loops = F)
    
    z=cbind(p,Secluded)
    
    a=rbind(a,z)
    print(j)
  }
  
  names(a)=c("Den","Secluded")
  
  a$Den=as.numeric(as.character(a$Den))
  a$Den=round(a$Den,5)
  
  
  df$Density=a$Den
  df$Secluded=a$Secluded
  
  isolated_comm=subset(df,df$Size>2&df$Secluded=='YES')
  isolated_comm$Status='Not Required'
  
  
  df1=subset(df,df$Size>2&df$Secluded=='NO')
  
  x=as.data.frame(isolated_comm$Comm_label)
  y=as.data.frame(df1$Comm_label)
  
  names(x)='Comm_label'
  names(y)='Comm_label'
  filter_id=rbind(x,y)
  
  
  t = quantile(df1$Size[df1$Size > 6], c(.25, .5))
  
  df1$Size_label = cut(
    df1$Size,
    c(0, 6, t, max(df1$Size)),
    labels = c('Small World', 'Small', 'Medium', 'Large')
  )
  
  df1$Density_label[df1$Size_label == "Small"] =
    cut(
      df1$Density[df1$Size_label == "Small"],
      c(0, quantile(df1$Density[df1$Size_label == "Small"], c(0.5)), 1),
      labels = c("Below", "Above"),
      include.lowest = F,
      right = T
    )
  
  
  df1$Density_label[df1$Size_label == "Medium"] =
    cut(
      df1$Density[df1$Size_label == "Medium"],
      c(0, quantile(df1$Density[df1$Size_label == "Medium"], c(0.5)), 1),
      labels = c("Below", "Above"),
      include.lowest = F,
      right = T
    )
  
  
  df1$Density_label[df1$Size_label == "Large"] =
    cut(
      df1$Density[df1$Size_label == "Large"],
      c(0, quantile(df1$Density[df1$Size_label == "Large"], c(0.5)), 1),
      labels = c("Below", "Above"),
      include.lowest = F,
      right = T
    )
  
  
  #Case 1, Small comm
  df1$Status[df1$Size_label == "Small" &
               df1$Density_label == "1"] = "Small & Loosely Knit"
  
  df1$Status[df1$Size_label == "Small" &
               df1$Density_label == "2"] = "Small & Closely Knit"
  
  #Case 2, Medium comm
  df1$Status[df1$Size_label == "Medium" &
               df1$Density_label == "1"] = "Medium & Loosely Knit"
  
  df1$Status[df1$Size_label == "Medium" &
               df1$Density_label == "2"] = "Medium & Closely Knit"
  
  #Case 3, Large comm
  df1$Status[df1$Size_label == "Large" &
               df1$Density_label == "1"] = "Large & Loosely Knit"
  
  df1$Status[df1$Size_label == "Large" &
               df1$Density_label == "2"] = "Large & Closely Knit"
  
  #Case 4, Small World
  df1$Status[df1$Size_label == "Small World"] ="Small World"
  
  
  G2_density=graph.density(G2,loops = F)
  
  df1$Status[df1$Density <= G2_density] = "VERY LOW DENSITY"
  
  
  
  nextlevel <- df1 %>% filter(df1$Status == "Large & Loosely Knit" |
                                df1$Status == "VERY LOW DENSITY")%>%select(Comm_label)
  
  maxsize <- max(df1$Size)
  
  density_thresholds <- df1 %>% 
    group_by(Size_label) %>% 
    summarise(thresh = quantile(Density, 0.5))
  
  df1=df1 %>% filter(!Comm_label  %in% nextlevel$Comm_label)
  
  nodes_df <- get.data.frame(G,"vertices") %>%
    select(name)%>%
    mutate(Comm_label=membership(Comm_obj)) %>%
    filter(Comm_label  %in% filter_id$Comm_label)
  
  nodes_df$Comm_label=as.character(nodes_df$Comm_label)
  
  
  sub_df=data.frame(
  )
  
  for (i in 1:length(nextlevel$Comm_label)) {
    
    v=as.numeric(nextlevel$Comm_label[i])
    C1=induced.subgraph(G, which(membership(Comm_obj)==as.numeric(nextlevel$Comm_label[i])))
    CD_2=cluster_louvain(C1,weights = E(C1)$weight)
    
    df2=as.data.frame(sizes(CD_2))
    names(df2)=c("Comm_label","Size")
    
    df2=subset(df2,df2$Size>2)
    
    b=data.frame(
      numeric()
    )
    
    
    for(k in 1:length(df2$Comm_label)){
      C2=induced.subgraph(C1, which(membership(CD_2)==as.numeric(as.character(df2$Comm_label[k]))))
      p = graph.density(C2,loops = F)
      b=rbind(b,p)
    }
    
    names(b)="Den"
    b$Den=round(b$Den,5)
    
    df2$Density=b$Den
    
    df2$Comm_label=paste0(v, sep="_", df2$Comm_label)
    
    df2$Secluded='NO'
    
    sub_df=bind_rows(sub_df,df2)
    
    
    
    nodes_df2 <- get.data.frame(C1,"vertices") %>%
      select(name)%>%
      mutate(Comm_label=membership(CD_2))%>%
      mutate(Comm_label = paste0(v, sep="_", Comm_label))%>%
      filter(Comm_label  %in% df2$Comm_label)
    
    nodes_df = bind_rows(nodes_df, nodes_df2)
    
    print(i)
    
  }
  
  
  
  sub_df$Size_label = cut(sub_df$Size,
                          c(0, 6, t, maxsize),
                          labels = c('Small World', 'Small', 'Medium', 'Large'))
  

  
  #sub_df$Status <- rep(NA, nrow(sub_df))
  sub_df$Status='a'
  
  sub_df$Status[sub_df$Size_label == "Small"] = ifelse(
    sub_df$Density[sub_df$Size_label == "Small"] < 
      density_thresholds$thresh[2],
    "Small & Loosely Knit",
    "Small & Closely Knit"
  )
  
  sub_df$Status[sub_df$Size_label == "Medium"] = ifelse(
    sub_df$Density[sub_df$Size_label == "Medium"] < 
      density_thresholds$thresh[3],
    "Medium & Loosely Knit",
    "Medium & Closely Knit"
  )
  
  sub_df$Status[sub_df$Size_label == "Large"] = ifelse(
    sub_df$Density[sub_df$Size_label == "Large"] < 
      density_thresholds$thresh[4],
    "Large & Loosely Knit",
    "Large & Closely Knit"
  )
  
  #Case 4, Small World
  sub_df$Status[sub_df$Size_label == "Small World"] = "Small World"
  
  
  sub_df=sub_df%>% select(Comm_label,Size,Density,Secluded,Status)
  
  df1=df1%>% select(Comm_label,Size,Density,Secluded,Status)
  df1$Comm_label=as.character(df1$Comm_label)
  df1$Secluded=as.character(df1$Secluded)
  
  sub_df$Comm_label=as.character(sub_df$Comm_label)
  
  isolated_comm$Secluded=as.character(isolated_comm$Secluded)
  isolated_comm$Comm_label=as.character(isolated_comm$Comm_label)
  
  
  Comm_details=bind_rows(isolated_comm,df1,sub_df)
  
  final_nodes=inner_join(nodes_df,Comm_details,by='Comm_label')
  
  #rithwik <<-Comm_details
  #vipin<<-final_nodes
  
  #write.csv(Comm_details,"/disk1/Workspace/SNA_maxis_oct/Sabah/Sabah_Comm_details.csv",row.names = F)
  #write.csv(Comm_details,"/Users/vipin.chauhan/Desktop/TMFDocs/Sabah_Comm_details.csv",row.names = F)
  write.csv(Comm_details,"/disk1/Workspace/SNA_maxis_oct/KL_Comm_details.csv",row.names = F)
  #write.csv(final_nodes,"/disk1/Workspace/SNA_maxis_oct/Sabah/Sabah_nodes_comm.csv",row.names = F)
  #write.csv(Comm_details,"/Users/vipin.chauhan/Desktop/TMFDocs/Sabah_nodes_comm.csv",row.names = F)
  write.csv(final_nodes,"/disk1/Workspace/SNA_maxis_oct/KL_nodes_comm.csv",row.names = F)
  
}

Node_Comm_fn(G2,Fast_community_unfolding)

