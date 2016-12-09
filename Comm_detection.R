library(dplyr)
library(igraph)
library(parallel)
options(scipen = 999)


#edges_df=edges_df %>% select(max_sna_opn_onnet.msisdn,max_sna_opn_onnet.other_party_number,max_sna_opn_onnet.total_calls,max_sna_opn_onnet.total_duration,max_sna_opn_onnet.total_sms,max_sna_opn_onnet.distinct_days)

names(edges_df)=c("CALLER_NUMBER","CALLED_NUMBER","TOTAL_CALLS","TOTAL_DURATION","TOTAL_SMS","DISTINCT_DAYS")
names(nodes_df)=c("msisdn","home_loc","is_post_paid","is_mop","is_churn","is_acquisition","nationality","arpu","package_desc","total_calls","total_outgoing_calls","total_incoming_calls","total_calls_roaming","total_duration","total_duration_roaming","total_outgoing_sms","total_flux","roaming_total_flux","night_calls","night_duration")

#Converting to character for measuring length
edges_df$CALLER_NUMBER=as.character(edges_df$CALLER_NUMBER)
edges_df$CALLED_NUMBER=as.character(edges_df$CALLED_NUMBER)

#edges_df=subset(edges_df,nchar(edges_df$CALLER_NUMBER,type = "char")>6&nchar(edges_df$CALLED_NUMBER,type = "char")>6)
#edges_df=subset(edges_df,nchar(edges_df$CALLER_NUMBER,type = "char")<15&nchar(edges_df$CALLED_NUMBER,type = "char")<15)

#aggregating data
start.time=Sys.time()
G1=simplify(graph.data.frame(edges_df,directed = F,vertices=nodes_df),remove.loops = TRUE,edge.attr.comb=list(TOTAL_CALLS="sum",TOTAL_DURATION="sum",TOTAL_SMS="sum",DISTINCT_DAYS="max")) 
rm(nodes_df)

end.time=Sys.time()
agg.data=round(difftime(end.time,start.time,units = "mins"),5)
cat( "agg.data=",agg.data,"\n")


rm(edges_df)


#Degree
V(G1)$Degree=degree(G1, v=V(G1),mode=c("all"),loops=F)

outlier_limit=quantile(V(G1)$Degree,probs=0.999)

#Remove outlier node with exceptionally high degree
a=V(G1)$name[degree(G1)>=outlier_limit[1]]
write.csv(a,file="/disk1/Workspace/SNA_maxis/Data/Outlier_Phno_Deg_KL.csv")
rm(a)

G1=delete.vertices(G1,c(V(G1)$name[degree(G1)>=outlier_limit[1]])) #WE CAN VARY THIS

V(G1)$Degree=degree(G1, v=V(G1),mode=c("all"),loops=F)

G1=delete.vertices(G1,which(degree(G1)==0))


#Filtered edgelist
detach(package:dplyr)
edges_df=as_data_frame(G1, what="edges")
nodes_df=as_data_frame(G1, what="vertices")

rm(G1)
edges_df=subset(edges_df,edges_df$TOTAL_DURATION<=2678400)
edges_df=subset(edges_df,edges_df$DISTINCT_DAYS<=31)

#Binnng
edges_df$TOTAL_CALLS=as.numeric(cut(edges_df$TOTAL_CALLS,breaks = c(-1,1,4,10,50,Inf),labels = c(1,2,3,4,5),right = T,include.lowest = F))
edges_df$TOTAL_DURATION=as.numeric(cut(edges_df$TOTAL_DURATION,breaks = c(-1,10,60,180,600,Inf),labels = c(1,2,3,4,5),right = T,include.lowest = F))
edges_df$TOTAL_SMS=as.numeric(cut(edges_df$TOTAL_SMS,breaks = c(-1,2,29,Inf),labels = c(1,2,3),right = T,include.lowest = F))
edges_df$TOTAL_SMS=log(edges_df$TOTAL_SMS)
edges_df$DISTINCT_DAYS=as.numeric(cut(edges_df$DISTINCT_DAYS,breaks = c(-1,3,6,20,Inf),labels = c(1,2,3,4),right = T,include.lowest = F))

#Scaling
edges_df$TOTAL_CALLS=round(1+(edges_df$TOTAL_CALLS - min(edges_df$TOTAL_CALLS))*9/diff(range(edges_df$TOTAL_CALLS)),5)
edges_df$TOTAL_DURATION=round(1+(edges_df$TOTAL_DURATION - min(edges_df$TOTAL_DURATION))*9/diff(range(edges_df$TOTAL_DURATION)),5)
edges_df$TOTAL_SMS=round(1+(edges_df$TOTAL_SMS - min(edges_df$TOTAL_SMS))*9/diff(range(edges_df$TOTAL_SMS)),5)
edges_df$DISTINCT_DAYS=round(1+(edges_df$DISTINCT_DAYS - min(edges_df$DISTINCT_DAYS))*9/diff(range(edges_df$DISTINCT_DAYS)),5)


#Calculating edge strength and again scaling (1-10)
edges_df$ES1=(edges_df$TOTAL_CALLS*0.25+edges_df$TOTAL_DURATION*0.25+edges_df$TOTAL_SMS*0.1+edges_df$DISTINCT_DAYS*0.4)

library(dplyr)
edges_df=edges_df %>% select(from,to,ES1)


detach(package:dplyr)
start.time=Sys.time()
#making weighted igraph object from refined edgelist 
G2=graph.data.frame(edges_df[,c(1,2)],directed = F,vertices=nodes_df)
E(G2)$weight=as.numeric(edges_df[,3]) 

	
end.time=Sys.time()
wt.graph.time=round(difftime(end.time,start.time,units = "mins"),5)
cat( "wt.graph.time=",wt.graph.time,"\n")


rm(edges_df)
V(G2)$Degree=degree(G2, v=V(G2),mode=c("all"),loops=F)

#Applying community detection algorithm
start.time=Sys.time()
Fast_community_unfolding=cluster_louvain(G2,weights = E(G2)$weight)

end.time=Sys.time()
CD_level1.data=round(difftime(end.time,start.time,units = "mins"),5)
cat( "CD_level1.data=",CD_level1.data,"\n")                                         

b=modularity(G2,membership(Fast_community_unfolding),weights = E(G2)$weight)#Greater the better (max 1)
cat( "modularity=",b,"\n")  







#Getting communtiy attributes for filtering
start.time=Sys.time()


df=as.data.frame(sizes(Fast_community_unfolding))

df=subset(df,df$Freq>2,select=c(1))

names(df)=c("Comm_label")


Community_df = data.frame(
Label=factor(),
Size=numeric(),
Pop_percent=numeric(),
Density=numeric(),
Secluded=character(),
MOP_perc=numeric(),
Post_paid_perc=numeric(),
Dominating_state=character(),
Dominating_plan=character(),
Plan_count=numeric(),
AVG_ARPU_Pre=numeric(),
AVG_ARPU_Post=numeric(),
Churn_perc=numeric(),
Aquired_perc=numeric(),
bangla_count=numeric(),
NPL_count=numeric(),
IND_count=numeric(),
INS_count=numeric(),
Outgng_calls_perc=numeric(),
Roaming_perc=numeric(),
Avg_Out_SMS=numeric(),
Avg_flux=numeric(),
Roaming_flux_perc=numeric(),
Night_dur_Perc=numeric()
)

for (i in as.numeric(df$Comm_label)) {
  
  C1=induced.subgraph(G2, which(membership(Fast_community_unfolding)==i))
  mean_deg_bfr=mean(V(C1)$Degree)
  V(C1)$Degree=degree(C1, v=V(C1),mode=c("all"),loops=F)
  Label=i
  Size=vcount(C1)
  Pop_percent=round(vcount(C1)*100/vcount(G2),3)
  Density=as.data.frame(graph.density(C1,loops = F))
  #Trans=as.data.frame(transitivity(C1,type = "global"))
  #Median_Degree=median(V(C1)$Degree)
  Secluded=ifelse((mean_deg_bfr-mean(V(C1)$Degree))==0,"YES","NO")
  
  MOP_cnt=sum(V(C1)$is_mop==1|V(C1)$is_mop==2,na.rm=T)
  Post_paid_cnt=sum(V(C1)$is_post_paid==1,na.rm=T)
  state=as.data.frame(table(V(C1)$home_loc))
  state=state[order(-state$Freq),]
  Dominating_state=state$Var1[1]

plan=as.data.frame(table(V(C1)$package_desc))

 if (nrow(plan)>0){
 plan=plan[order(-plan$Freq),]
 Dominating_plan=as.character(plan$Var1[1])
 Plan_count=plan$Freq[1]
}else {
   Dominating_plan="No Plan"
  Plan_count=0
 }
TTL_ARPU_Pre=sum(V(C1)$arpu[V(C1)$is_post_paid==0],na.rm=T)
TTL_ARPU_MOP=sum(V(C1)$arpu[V(C1)$is_mop==1|V(C1)$is_mop==2],na.rm=T)
TTL_APRU_non_MOP=sum(V(C1)$arpu[V(C1)$is_mop==0|V(C1)$is_post_paid==1],na.rm=T)

 # AVG_ARPU=round(sum(V(C1)$arpu,na.rm=T)/(vcount(C1)-sum(is.na(V(C1)$arpu))),3)
#  Churn_perc=round(sum(V(C1)$is_churn==1,na.rm=T)*100/Size,3)
 # Aquired_perc=round(sum(V(C1)$is_acquisition==1,na.rm=T)*100/Size,3)
  bangla_count=sum(V(C1)$nationality=='BNG',na.rm=T)
  NPL_count=sum(V(C1)$nationality=='NPL',na.rm=T)
  IND_count=sum(V(C1)$nationality=='IND',na.rm=T)
  INS_count=sum(V(C1)$nationality=='INS',na.rm=T)
  Outgng_calls_perc=round(sum(V(C1)$total_outgoing_calls,na.rm=T)*100/sum(V(C1)$total_calls,na.rm=T),3)
  Roaming_perc=round(sum(V(C1)$total_calls_roaming,na.rm=T)*100/sum(V(C1)$total_calls,na.rm=T),3)
  Avg_Out_SMS=round(sum(V(C1)$total_outgoing_sms,na.rm=T)/(vcount(C1)-sum(is.na(V(C1)$total_outgoing_sms))),3)
  Avg_flux=round(sum(V(C1)$total_flux,na.rm=T)/(vcount(C1)-sum(is.na(V(C1)$total_flux))),1)
  Roaming_flux_perc=round(sum(V(C1)$roaming_total_flux,na.rm=T)*100/sum(V(C1)$total_flux,na.rm=T),3)
  Night_dur_Perc=round(sum(V(C1)$night_duration,na.rm=T)*100/sum(V(C1)$total_duration,na.rm=T),3)

 a=cbind(Label,Size,Pop_percent,Density,Secluded,MOP_cnt,Post_paid_cnt,Dominating_state,Dominating_plan,Plan_count,TTL_ARPU_Pre,TTL_ARPU_MOP,TTL_APRU_non_MOP,bangla_count,NPL_count,IND_count,INS_count,Outgng_calls_perc,Roaming_perc,Avg_Out_SMS,Avg_flux,Roaming_flux_perc,Night_dur_Perc)
  
Community_df = rbind(Community_df, a)

  
 comm_attr_perc=nrow(Community_df)*100/nrow(df)
  cat( "comm_attr_perc_covered=",comm_attr_perc,"\n")
   
}
  

  rm(a,Label,Size,Pop_percent,Density,i,C1,mean_deg_bfr,Secluded,state,Dominating_state,Plan_count,Post_paid_cnt,MOP_cnt,Dominating_plan,TTL_ARPU_Pre,TTL_ARPU_MOP,TTL_APRU_non_MOP,bangla_count,NPL_count,IND_count,INS_count,Outgng_calls_perc,Roaming_perc,Avg_Out_SMS,Avg_flux,Roaming_flux_perc,Night_dur_Perc)


rm(df)

names(Community_df)=c("Community Label","Size","Pop_percent","Density","Secluded","MOP_cnt","Post_paid_cnt","Dominating_state","Dominating_plan","Plan_count","TTL_ARPU_Pre","TTL_ARPU_MOP","TTL_APRU_non_MOP","bangla_count","NPL_count","IND_count","INS_count","Outgng_calls_perc","Roaming_perc","Avg_Out_SMS","Avg_flux","Roaming_flux_perc","Night_dur_Perc")

    Community_df$`Community Label`=as.factor(Community_df$`Community Label`)


#Rounding off
Community_df$Density=round(Community_df$Density,5)
#Community_df$Transitivity=round(Community_df$Transitivity,2)
Community_df=Community_df[order(-Community_df$Size),]


#Separating isolated & non isolated communities
isolated_comm=subset(Community_df,Community_df$Secluded=="YES")
Non_isolated_comm=subset(Community_df,Community_df$Secluded=="NO")
rm(Community_df)


#Breaking size and density on quantiles

t = quantile(Non_isolated_comm$Size[Non_isolated_comm$Size > 6], c(.25, .5))

Non_isolated_comm$Size_label=cut(Non_isolated_comm$Size, c(0, 6, t, max(Non_isolated_comm$Size)), labels = c('Small World', 'Small', 'Medium', 'Large'))



Non_isolated_comm$Density_label[Non_isolated_comm$Size_label=="Small"]=cut(Non_isolated_comm$Density[Non_isolated_comm$Size_label=="Small"],c(0,quantile(Non_isolated_comm$Density[Non_isolated_comm$Size_label=="Small"],c(0.5)),1),labels = c("Below","Above"),include.lowest = F,right=T)

Non_isolated_comm$Density_label[Non_isolated_comm$Size_label=="Medium"]=cut(Non_isolated_comm$Density[Non_isolated_comm$Size_label=="Medium"],c(0,quantile(Non_isolated_comm$Density[Non_isolated_comm$Size_label=="Medium"],c(0.5)),1),labels = c("Below","Above"),include.lowest = F,right=T)

Non_isolated_comm$Density_label[Non_isolated_comm$Size_label=="Large"]=cut(Non_isolated_comm$Density[Non_isolated_comm$Size_label=="Large"],c(0,quantile(Non_isolated_comm$Density[Non_isolated_comm$Size_label=="Large"],c(0.5)),1),labels = c("Below","Above"),include.lowest = F,right=T)
    




#Naming each Non_isolated community

#Case 1, Small comm
Non_isolated_comm$Status[Non_isolated_comm$Size_label=="Small"&
                           Non_isolated_comm$Density_label=="1"]="Small & Loosely Knit"

Non_isolated_comm$Status[Non_isolated_comm$Size_label=="Small"&
                           Non_isolated_comm$Density_label=="2"]="Small & Closely Knit"

#Case 2, Medium comm
Non_isolated_comm$Status[Non_isolated_comm$Size_label=="Medium"&
                           Non_isolated_comm$Density_label=="1"]="Medium & Loosely Knit"

Non_isolated_comm$Status[Non_isolated_comm$Size_label=="Medium"&
                           Non_isolated_comm$Density_label=="2"]="Medium & Closely Knit"

#Case 3, Large comm
Non_isolated_comm$Status[Non_isolated_comm$Size_label=="Large"&
                           Non_isolated_comm$Density_label=="1"]="Large & Loosely Knit"

Non_isolated_comm$Status[Non_isolated_comm$Size_label=="Large"&
                           Non_isolated_comm$Density_label=="2"]="Large & Closely Knit"

#Case 4, Small World
Non_isolated_comm$Status[Non_isolated_comm$Size_label=="Small World"]="Small World"


#Case 5, communties with density lesser than overall graph
G2_density=graph.density(G2,loops = F)
Non_isolated_comm$Status[Non_isolated_comm$Density<=G2_density]="VERY LOW DENSITY"
rm(G2_density)



end.time=Sys.time()
comm.attr.data=round(difftime(end.time,start.time,units = "mins"),5)
cat("comm.attr.data=",comm.attr.data,"\n")






start.time=Sys.time()
Cycle1=subset(Non_isolated_comm,Non_isolated_comm$Status=="Large & Loosely Knit"|Non_isolated_comm$Status=="VERY LOW DENSITY",select = c(1))


#Again dividing select communties to increase density
subcomm_df = data.frame(
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

for (i in Cycle1$`Community Label`) {
  Original_Comm=i
  C2=induced.subgraph(G2, which(membership(Fast_community_unfolding)==i))
  CD_2=cluster_louvain(C2,weights = E(C2)$weight)
  #Sub_comm=unique(CD_2$membership)

df=as.data.frame(sizes(CD_2))

df=subset(df,df$Freq>2,select=c(1))

names(df)=c("Comm_label")


  for (j in df$Comm_label) {
    C3=induced.subgraph(C2, which(membership(CD_2)==j))
    Size=vcount(C3)
    Pop_percent=round(vcount(C3)*100/vcount(G2),3)
    Density=as.data.frame(graph.density(C3,loops = F))
    #Trans=as.data.frame(transitivity(C3,type = "global"))
    #Median_Degree=median(V(C3)$Degree)
    #Value=sum(E(C3)$TOTAL_VALUE)
  
  MOP_cnt=sum(V(C3)$is_mop==1|V(C3)$is_mop==2,na.rm=T)
  Post_paid_cnt=sum(V(C3)$is_post_paid==1,na.rm=T)
  state=as.data.frame(table(V(C3)$home_loc))
  state=state[order(-state$Freq),]
  Dominating_state=state$Var1[1]


 plan=as.data.frame(table(V(C3)$package_desc))

   if (nrow(plan)>0){
         plan=plan[order(-plan$Freq),]
           Dominating_plan=plan$Var1[1]
             Plan_count=plan$Freq[1]
   }else {
           Dominating_plan=NA
               Plan_count=0
                   }

TTL_ARPU_Pre=sum(V(C1)$arpu[V(C1)$is_post_paid==0],na.rm=T)
TTL_ARPU_MOP=sum(V(C1)$arpu[V(C1)$is_mop==1|V(C1)$is_mop==2],na.rm=T)
TTL_APRU_non_MOP=sum(V(C1)$arpu[V(C1)$is_mop==0|V(C1)$is_post_paid==1],na.rm=T)


  #AVG_ARPU=round(sum(V(C3)$arpu,na.rm=T)/(vcount(C3)-sum(is.na(V(C3)$arpu))),3)
 # Churn_perc=round(sum(V(C3)$is_churn==1,na.rm=T)*100/Size,3)
 # Aquired_perc=round(sum(V(C3)$is_acquisition==1,na.rm=T)*100/Size,3)
  bangla_count=sum(V(C3)$nationality=='BNG',na.rm=T)
  NPL_count=sum(V(C3)$nationality=='NPL',na.rm=T)
  IND_count=sum(V(C3)$nationality=='IND',na.rm=T)
  INS_count=sum(V(C3)$nationality=='INS',na.rm=T)
  Outgng_calls_perc=round(sum(V(C3)$total_outgoing_calls,na.rm=T)*100/sum(V(C3)$total_calls,na.rm=T),3)
  Roaming_perc=round(sum(V(C3)$total_calls_roaming,na.rm=T)*100/sum(V(C3)$total_calls,na.rm=T),3)
  Avg_Out_SMS=round(sum(V(C3)$total_outgoing_sms,na.rm=T)/(vcount(C3)-sum(is.na(V(C3)$total_outgoing_sms))),3)
  Avg_flux=round(sum(V(C3)$total_flux,na.rm=T)/(vcount(C3)-sum(is.na(V(C3)$total_flux))),3)
  Roaming_flux_perc=round(sum(V(C3)$roaming_total_flux,na.rm=T)*100/sum(V(C3)$total_flux,na.rm=T),3)
  Night_dur_Perc=round(sum(V(C3)$night_duration,na.rm=T)*100/sum(V(C3)$total_duration,na.rm=T),3)


    b=cbind(Original_Comm,j,Size,Pop_percent,Density,MOP_cnt,Post_paid_cnt,Dominating_state,Dominating_plan,Plan_count,TTL_ARPU_Pre,TTL_ARPU_MOP,TTL_APRU_non_MOP,bangla_count,NPL_count,IND_count,INS_count,Outgng_calls_perc,Roaming_perc,Avg_Out_SMS,Avg_flux,Roaming_flux_perc,Night_dur_Perc)
    subcomm_df = rbind(subcomm_df, b)
    rm(C3,j,Size,Pop_percent,Density,MOP_cnt,Post_paid_cnt,Dominating_state,Dominating_plan,Plan_count,TTL_ARPU_Pre,TTL_ARPU_MOP,TTL_APRU_non_MOP,bangla_count,NPL_count,IND_count,INS_count,Outgng_calls_perc,Roaming_perc,Avg_Out_SMS,Avg_flux,Roaming_flux_perc,Night_dur_Perc)
  }
 subcomm_attr_perc=length(unique(subcomm_df$Original_Comm))*100/nrow(Cycle1)
  cat( "sub_comm_attr_perc_covered=",subcomm_attr_perc,"\n")
  rm(subcomm_attr_perc)

  rm(C2,i,CD_2,b,Original_Comm)
}

names(subcomm_df)=c("Original_Community","Sub_Commnity","Size","Pop_percent","Density","MOP_cnt","Post_paid_cnt","Dominating_state","Dominating_plan","Plan_count","TTL_ARPU_Pre","TTL_ARPU_MOP","TTL_APRU_non_MOP","bangla_count","NPL_count","IND_count","INS_count","Outgng_calls_perc","Roaming_perc","Avg_Out_SMS","Avg_flux","Roaming_flux_perc","Night_dur_Perc")


#subcomm_df=merge(x = subcomm_df, y = Non_isolated_comm, by.x = "Original_Community", by.y="Community Label",all.x = T)

#library(dplyr)
#subcomm_df=subcomm_df %>% select(Original_Community,Sub_Commnity,Size_subcomm,Pop_percent_subcomm,Density_subcomm,Transitivity_subcomm,Median_Degree_subcomm,MOP_perc_sub,Post_paid_perc_sub,Dominating_state_sub,Density,)

#detach(package:dplyr)

#names(subcomm_df)=c("Original_Community","Sub-Commnity","Size","Pop_percent","Density","Transitivity","Median_Degree","MOP_perc","Post_paid_perc","Dominating_state","Density_Initially")


subcomm_df$Density=round(subcomm_df$Density,5)
#subcomm_df$Transitivity=round(subcomm_df$Transitivity,2)


#Breaking size and density on quantiles
subcomm_df$Size_label=cut(subcomm_df$Size, c(0, 6, t, max(Non_isolated_comm$Size)), labels = c('Small World', 'Small', 'Medium', 'Large'))

subcomm_df$Density_label[subcomm_df$Size_label=="Small"]=cut(subcomm_df$Density[subcomm_df$Size_label=="Small"],c(0,quantile(Non_isolated_comm$Density[Non_isolated_comm$Size_label=="Small"],c(0.5)),1),labels = c("Below","Above"),include.lowest = F,right=T)

subcomm_df$Density_label[subcomm_df$Size_label=="Medium"]=cut(subcomm_df$Density[subcomm_df$Size_label=="Medium"],c(0,quantile(Non_isolated_comm$Density[Non_isolated_comm$Size_label=="Medium"],c(0.5)),1),labels = c("Below Average","Above Average"),include.lowest = F,right=T)

subcomm_df$Density_label[subcomm_df$Size_label=="Large"]=cut(subcomm_df$Density[subcomm_df$Size_label=="Large"],c(0,quantile(Non_isolated_comm$Density[Non_isolated_comm$Size_label=="Large"],c(0.5)),1),labels = c("Below Average","Above Average"),include.lowest = F,right=T)



subcomm_df$Status[subcomm_df$Size_label=="Small"&
                           subcomm_df$Density_label=="1"]="Small & Loosely Knit"

subcomm_df$Status[subcomm_df$Size_label=="Small"&
                           subcomm_df$Density_label=="2"]="Small & Closely Knit"

#Case 2, Medium comm
subcomm_df$Status[subcomm_df$Size_label=="Medium"&
                           subcomm_df$Density_label=="1"]="Medium & Loosely Knit"

subcomm_df$Status[subcomm_df$Size_label=="Medium"&
                           subcomm_df$Density_label=="2"]="Medium & Closely Knit"

#Case 3, Large comm
subcomm_df$Status[subcomm_df$Size_label=="Large"&
                           subcomm_df$Density_label=="1"]="Large & Loosely Knit"

subcomm_df$Status[subcomm_df$Size_label=="Large"&
                           subcomm_df$Density_label=="2"]="Large & Closely Knit"

#Case 4, Small World
subcomm_df$Status[subcomm_df$Size_label=="Small World"]="Small World"

#Case 5
#subcomm_df$Status[subcomm_df$Size_label==<NA>]="Bad Split"

#subcomm_df$Size_label=cut(subcomm_df$Size,quantile(subcomm_df$Size,probs=c(0,0.5,1)),labels = c("Below Half","Above Half"),include.lowest = T)
#subcomm_df$Density_label=cut(subcomm_df$Density,quantile(subcomm_df$Density,probs=c(0,0.25,1)),labels = c("Low","OK"),include.lowest = T)



#Case 1, Small comm
#subcomm_df$Status[subcomm_df$Size_label=="Below Half"&
 #                   subcomm_df$Density_label=="Low"]="Small & Loosely Knit"

#subcomm_df$Status[subcomm_df$Size_label=="Below Half"&
 #                  subcomm_df$Density_label=="OK"]="Small & Closely Knit"

#Case 2, Bigger comm
#subcomm_df$Status[subcomm_df$Size_label=="Above Half"&
 #                   subcomm_df$Density_label=="Low"]="Decent Sized & Loosely Knit"

#subcomm_df$Status[subcomm_df$Size_label=="Above Half"&
 #                   subcomm_df$Density_label=="OK"]="Decent Sized & Closely Knit"


end.time=Sys.time()
subcomm.data=round(difftime(end.time,start.time,units = "mins"),5)
cat("subcomm.data=",subcomm.data,"\n")



save(Fast_community_unfolding,isolated_comm,Non_isolated_comm,subcomm_df,G2,file = "/disk1/Workspace/SNA_maxis_oct/CD_L1_KL.RData")


table(subcomm_df$Status)




	    
