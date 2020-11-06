#Develop correlation network of adverse events
setwd("........")
library(readr)
library(stringr)

nerve_out = read.delim("", header=FALSE)
nerve_df=read_csv("")

psychiatric_out = read.delim("", header=FALSE)
psychiatric_df=read_csv("")

tissue_out = read.delim("", header=FALSE)
tissue_df=read_csv("")

gas_out = read.delim("", header=FALSE)
gas_df=read_csv("")

#Get the frequency of pt level
#Top ten: Dizziness, Syncope, Headache, Loss of consciousness, Hypoaesthesia,Convulsion,Paraesthesia,Tremor,
#Muscular weakness, Dyskinesia
nerve_sum = as.data.frame(table(nerve_out$V1))
attach(nerve_sum)
nerve_sum2 = nerve_sum[order(-Freq),]
detach(nerve_sum)


pysch_sum = as.data.frame(table(psychiatric_out$V1))
attach(pysch_sum)
pysch_sum2 = pysch_sum[order(-Freq),]
detach(pysch_sum)

tissue_sum = as.data.frame(table(tissue_out$V1))
attach(tissue_sum)
tissue_sum2 = tissue_sum[order(-Freq),]
detach(tissue_sum)

gas_sum = as.data.frame(table(gas_out$V1))
attach(gas_sum)
gas_sum2 = gas_sum[order(-Freq),]
detach(gas_sum)


#List of Unique Nervous System Disorders Pt
nerve_pt = as.data.frame(nerve_sum2[,1])
nerve_df2 = nerve_df

pysch_pt = as.data.frame(pysch_sum2[,1])
psychiatric_df2 = psychiatric_df

tissue_pt = as.data.frame(tissue_sum2[,1])
tissue_df2 = tissue_df

gas_pt = as.data.frame(gas_sum2[,1])
gas_df2 = gas_df


#Create Pt-level wide format dummy variables

for (i in 1:nrow(nerve_pt)){
  nerve_df2[,i+2]=str_count(nerve_df2$Symptoms, as.character(nerve_pt[i,]))
  names(nerve_df2)[i+2]=as.character(nerve_pt[i,])
}



for (i in 1:nrow(pysch_pt)){
  psychiatric_df2[,i+2]=str_count(psychiatric_df2$Symptoms, as.character(pysch_pt[i,]))
  names(psychiatric_df2)[i+2]=as.character(pysch_pt[i,])
}

for (i in 1:nrow(tissue_pt)){
  tissue_df2[,i+2]=str_count(tissue_df2$Symptoms, as.character(tissue_pt[i,]))
  names(tissue_df2)[i+2]=as.character(tissue_pt[i,])
}

for (i in 1:nrow(gas_pt)){
  gas_df2[,i+2]=str_count(gas_df2$Symptoms, as.character(gas_pt[i,]))
  names(gas_df2)[i+2]=as.character(gas_pt[i,])
}


 


#Reduce multiple same pt/patient to one, need to ask why same pts occured
nerve_df3 = nerve_df2[,-2]

psychiatric_df3 = psychiatric_df2[,-2]

tissue_df3 = tissue_df2[,-2]

gas_df3 = gas_df2[,-2]

for (i in c(3:52)){
  nerve_df3[,i]=ifelse(nerve_df3[,i]>=1,1,0)
}

for (i in c(3:52)){
  psychiatric_df3[,i]=ifelse(psychiatric_df3[,i]>=1,1,0)
}


for (i in c(3:52)){
  tissue_df3[,i]=ifelse(tissue_df3[,i]>=1,1,0)
}


for (i in c(3:52)){
  gas_df3[,i]=ifelse(gas_df3[,i]>=1,1,0)
}

names(nerve_df3)[1] = "ID"


#Get Nerve-psychi Nerve-tissue Nerve-gas pairwise PT-term df
nerve_pyschi = merge(psychiatric_df3, nerve_df3, all= T,by ="ID",sort = T)

nerve_tissue = merge(tissue_df3, nerve_df3,all= T, by ="ID",sort = T)

nerve_gas = merge(gas_df3, nerve_df3, all= T,by ="ID", sort = T)

#Get the PT-term matrix only
nerve_pyschi2 = nerve_pyschi[,-1]
nerve_tissue2 = nerve_tissue[,-1]
nerve_gas2 = nerve_gas[,-1]

nerve_pyschi2[is.na(nerve_pyschi2)] = 0
nerve_tissue2[is.na(nerve_tissue2)] = 0
nerve_gas2[is.na(nerve_gas2)] = 0


#network
library(tidyverse)
library(corrr)
library(igraph)
library(ggraph)
library(psych)

#Get Correlation Matrix p-values, adjust for multiple testing
cor_nervegas = corr.test(nerve_gas2[,1:ncol(nerve_gas2)],nerve_gas2[,1:ncol(nerve_gas2)], adjust = "holm")
cor_nervegas_p = as.data.frame(cor_nervegas$p)
library(data.table)
setDT(cor_nervegas_p, keep.rownames = TRUE)[]
library(reshape2)
cor_nervegas_p2 = melt(cor_nervegas_p)
cor_nervegas_p3 = cor_nervegas_p2[,c(2,1,3)]


#library(lsr)
tidy_cors = nerve_gas2 %>% 
  correlate() %>% 
  stretch()

tidy_cors$p = cor_nervegas_p3$value
# Convert correlations stronger than some value
# to an undirected graph object
graph_cors = tidy_cors %>% 
  filter(abs(r) > 0.4 & p<0.05) %>% 
  graph_from_data_frame(directed = FALSE)

# Plot
netplot = ggraph(graph_cors)
netplot$data$SOC_type = ifelse(netplot$data$name %in% gas_pt[,1],"Gas","nerve")

tiff('Nerve_Gas_network.tiff', units="in", width=12, height=9, res=300, compression = 'lzw')
netplot +geom_edge_link(aes(edge_alpha = abs(r), edge_width = abs(r), color = r)) +
  guides(edge_alpha = "none", edge_width = "none") +
  scale_edge_colour_gradientn(limits = c(-1, 1), colors = c("firebrick2", "purple","dodgerblue2")) +
  geom_node_point(aes(color = "black", size = 1, shape = SOC_type)) +
  geom_node_text(aes(label = name), repel = TRUE, size =2) +
  theme_graph() +
  labs(title = "Correlation Network of nerve_gas2 PT Terms with Moderate or above (r>0.4) correlation strength,p<0.05 after adjustment")

dev.off()

#Nerve Tissue Network

cor_nervetis = corr.test(nerve_tissue2[,1:ncol(nerve_tissue2)],nerve_tissue2[,1:ncol(nerve_tissue2)], adjust = "holm")
cor_nervetis_p = as.data.frame(cor_nervetis$p)
library(data.table)
setDT(cor_nervetis_p, keep.rownames = TRUE)[]
library(reshape2)
cor_nervetis_p2 = melt(cor_nervetis_p)
cor_nervetis_p3 = cor_nervetis_p2[,c(2,1,3)]


#library(lsr)
tidy_cors = nerve_tissue2 %>% 
  correlate() %>% 
  stretch()

tidy_cors$p = cor_nervetis_p3$value
# Convert correlations stronger than some value
# to an undirected graph object
graph_cors = tidy_cors %>% 
  filter(abs(r) > 0.2 & p<0.05) %>% 
  graph_from_data_frame(directed = FALSE)

# Plot
netplot = ggraph(graph_cors)
#Group by PT-terms, create new indicator variable SOC_Type
netplot$data$SOC_type = ifelse(netplot$data$name %in% tissue_pt[,1],"Tissue","nerve")

tiff('Nerve_tissue_network2.tiff', units="in", width=12, height=9, res=300, compression = 'lzw')
netplot +geom_edge_link(aes(edge_alpha = abs(r), edge_width = abs(r), color = r)) +
  guides(edge_alpha = "none", edge_width = "none") +
  scale_edge_colour_gradientn(limits = c(-1, 1), colors = c("firebrick2", "purple","dodgerblue2")) +
  geom_node_point(aes(color = "black", size = 1, shape = SOC_type)) +
  geom_node_text(aes(label = name), repel = TRUE, size =2) +
  theme_graph() +
  labs(title = "Correlation Network of nerve_tissue2 PT Terms with Moderate or above (r>0.2) correlation strength,p<0.05 after adjustment")

dev.off()

test_df =tidy_cors %>% 
  filter(abs(r) > 0.2 & p<0.05)



#Nerve pyschi
cor_nervepychi = corr.test(nerve_pyschi2[,1:ncol(nerve_pyschi2)],nerve_pyschi2[,1:ncol(nerve_pyschi2)], adjust = "holm")
cor_nervepychi_p = as.data.frame(cor_nervepychi$p)
library(data.table)
setDT(cor_nervepychi_p, keep.rownames = TRUE)[]
library(reshape2)
cor_nervepychi_p2 = melt(cor_nervepychi_p)
cor_nervepychi_p3 = cor_nervepychi_p2[,c(2,1,3)]


#library(lsr)
tidy_cors = nerve_pyschi2 %>% 
  correlate() %>% 
  stretch()

tidy_cors$p = cor_nervepychi_p3$value
# Convert correlations stronger than some value
# to an undirected graph object
graph_cors = tidy_cors %>% 
  filter(abs(r) > 0.2 & p<0.05) %>% 
  graph_from_data_frame(directed = FALSE)

# Plot
netplot = ggraph(graph_cors)
#Group by PT-terms, create new indicator variable SOC_Type
netplot$data$SOC_type = ifelse(netplot$data$name %in% tissue_pt[,1],"Tissue","nerve")

tiff('Nerve_tissue_network2.tiff', units="in", width=12, height=9, res=300, compression = 'lzw')
netplot +geom_edge_link(aes(edge_alpha = abs(r), edge_width = abs(r), color = r)) +
  guides(edge_alpha = "none", edge_width = "none") +
  scale_edge_colour_gradientn(limits = c(-1, 1), colors = c("firebrick2", "purple","dodgerblue2")) +
  geom_node_point(aes(color = "black", size = 1, shape = SOC_type)) +
  geom_node_text(aes(label = name), repel = TRUE, size =2) +
  theme_graph() +
  labs(title = "Correlation Network of nerve_tissue2 PT Terms with Moderate or above (r>0.2) correlation strength,p<0.05 after adjustment")

dev.off()

test_df =tidy_cors %>% 
  filter(abs(r) > 0.2 & p<0.05)





























##############################SOC Level#########################################
#chisq test
tbl = table(nerve$`Loss of consciousness`, nerve$SEX) 
tbl
chisq.test(tbl) 

tbl = table(nerve$Hypoaesthesia, nerve$SEX) 
tbl
chisq.test(tbl)

setwd("Y:/HVP_Project")
library(readr)
library(stringr)
#Clustering of SOCs level


library(readr)
library(stringr)
SOC_text = read_csv("ALLresultHPV918.csv")

HPV_ALL = read_csv("HPV_ALL.csv")


df = merge(SOC_text, HPV_ALL, by ="VAERS_ID", all.x =T, sort = T)
df2 = df[,c(match("VAERS_ID", names(df)),match("SOC", names(df)),
            match("AGE_YRS", names(df)),match("SEX", names(df))
)]

#List of SOC
SOC_list = c(
  "infections and infestations",
  "Neoplasms",
  "blood and lymphatic system disorders",
  "immune system disorders",
  "endocrine disorders",
  "metabolism and nutrition disorders",
  "psychiatric disorders",
  "nervous system disorders",
  "eye disorders",
  "ear and labyrinth disorders",
  "cardiac disorders",
  "vascular disorders",
  "respiratory",
  "gastrointestinal disorders",
  "hepatobiliary disorders",
  "skin and subcutaneous tissue disorders",
  "musculoskeletal and connective tissue disorders",
  "renal and urinary disorders",
  "pregnancy",
  "reproductive system and breast disorders",
  "congenital_familial and genetic disorders",
  "general disorders and administration site conditions",
  "investigations",
  "Injury_poisoning and procedural complications",
  "surgical and medical procedures",
  "social circumstances"
)


#Map and create new variables
soclist = list()
for (i in 1:26){
  socname = paste("SOC",i,sep="")
  soclist[i] = socname
  df2[,i+4] = str_count(df2$SOC, SOC_list[i])
  
  names(df2)[i+4]=soclist[i]
}

df2_soc_dummy = df2

df2$ADE = apply(df2[,c(match("SOC1",names(df2)):match("SOC26",names(df2)))],1,sum)

#Delete those without any soc
df3 = df2[df2$ADE>0,]
df4 = df3[!df3$SEX=="U",]


for (i in c(5:30)){
  names(df4)[i]=SOC_list[i-4]
}

colnames(df4)[23]="pregnancy, puerperium and perinatal conditions"
SOC_matrix = df4[,c(5:30)]
SOC_corr = cor(SOC_matrix)

dissimilarity = dist(SOC_corr)
dissimilarity = 1 - abs(SOC_corr)

distance = as.dist(dissimilarity)
clust = plot(hclust(distance), 
             main ="",xlab="", sub = "")



#Determine the optimal number of clustering
#Average silhouette method
library(factoextra)
require(cluster)


fviz_nbclust(SOC_corr, hcut, method = "silhouette", diss = distance,
             hc_method = "average") +
  geom_vline(xintercept = 3, linetype = 2)

clust= hclust(distance) 
plot(clust,main ="",xlab="", sub = "")

clust
rect.hclust(clust, k = 3, border = 2:4)



clust = hclust(distance, method = "average")
clust_cut = as.data.frame(cutree(clust, k = 3))

library(data.table)
setDT(clust_cut,keep.rownames = TRUE)[]
clust_cut = clust_cut[order(clust_cut$`cutree(clust, k = 3)`),]
clust1 = as.data.frame(clust_cut[clust_cut$`cutree(clust, k = 3)`==1,1])
clust2 = as.data.frame(clust_cut[clust_cut$`cutree(clust, k = 3)`==2,1])
clust3 = as.data.frame(clust_cut[clust_cut$`cutree(clust, k = 3)`==3,1])

clust1 = clust1[,1]
clust2 = clust2[,1]
clust3 = clust3[,1]





#Correlation Network of SOCs
library(tidyverse)
library(corrr)
library(igraph)
library(ggraph)
tidy_cors_SOC = SOC_matrix %>% 
  correlate() %>% 
  stretch()



graph_cors_SOC= tidy_cors_SOC %>% 
  filter() %>% 
  graph_from_data_frame(directed = FALSE)



netplot = ggraph(graph_cors_SOC)

netplot$data$clustN = ifelse(netplot$data$name %in% clust1,"Cluster 1",
                             ifelse(netplot$data$name %in% clust2,"Cluster2","Cluster3"))

tiff("SOCs_network.tiff", units = 'in', width = 12, height = 9, res = 300,
     compression = "lzw")
netplot +geom_edge_link(aes(edge_alpha = abs(r), edge_width = abs(r), color = r)) +
  guides(edge_alpha = "none", edge_width = "none") +
  scale_edge_colour_gradientn(limits = c(-1, 1), colors = c("firebrick2", "white","dodgerblue2")) +
  geom_node_point(aes(color = clustN, size = 4)) +
  geom_node_text(aes(label = name), repel = TRUE, size =5) +
  theme_graph() +
  labs(title ="")

dev.off()






#Get the correlations with nervous disorders
for (i in c(2:26)){
  cor_soc = cor(test[,1],test[,i])
  print(paste("correlation between", names(test)[1],"and", names(test)[i], "=", round(cor_soc,3)))
}

