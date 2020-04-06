#####################################
#####
#####
#####MIRNA DATA ANALYSIS
#####By: Darryl Nousome
#####
#####
#####################################

library(tidyverse)
library(readxl)
library(DESeq2)
library(sva)

#library(visNetwork)


select<-dplyr::select
#####Load in all of the data
##IDs first
setwd("~/shares/F/Prostate/miRNA/")



###Read in miRNA
##Batch 1:3::CALL BATCH 1
d<-list.files("Data/data_1",full.names = T)
mirnas_b1<-lapply(d,function(x){
  temp<-read_xls(x,skip=8,sheet=2)
  temp<-temp[-1:-2,]
  names(temp)<-gsub("[_-]",".",names(temp))
  names(temp)<-gsub("run101.","",names(temp))
  
  #rownames(temp)<-unlist(temp[,1])
  temp[,-1]
})

##Remove .1 in names of B2
names(mirnas_b1[[2]])<-sapply(strsplit(names(mirnas_b1[[2]]),"[.]"),function(x)paste(x[-length(x)],collapse="."))






###Batch 4:7 
b4<-read_xls("Data/data_2/100518_HTG_Ram-Prostate-FFPE_1.xls",skip=10,sheet=2) %>% select(2:24)
b5<-read_xls("Data/data_2/100818_HTG_Ram-Prostate-FFPE_2.xls",skip=10,sheet=2) %>% select(2:24) 
b6<-read_xls("Data/data_2/110518_HTG_Ram-Prostate-FFPE-3.xls",skip=10,sheet=2) %>% select(2:24) 
b7<-read_xls("Data/data_2/030619_HTG_Ram-Prostate-FFPE-4.xls",skip=10,sheet=2) %>% select(2:25) 

mirnas.b2<-cbind(b4,b5,b6,b7)




###Batch 8:9 Normal
normal_p<-read_xls("Data/data_3/159-071219_HTG_Ram-Prostate-FFPE-Normals-MiSeq.xls",skip=10,sheet=1) %>% select(2:24) 
nospread_p<-read_xls("Data/data_3/172-111919_HTG-miRNA-Wilson-Prostate-Primaries-Not-Spread_n=21_Norm-MiSeq.xls",skip=10,sheet=1) %>% select(2:22) 

normal_p_ids_sheet=unlist(read_xls("Data/data_3/159-071219_HTG_Ram-Prostate-FFPE-Normals-MiSeq.xls")[8,-1])
normal_p_ids_sheet=sapply(strsplit(normal_p_ids_sheet,"[-_]"),function(x)paste0(x[2],"-",x[3]))


nospread_p_ids_sheet=unlist(read_xls("Data/data_3/172-111919_HTG-miRNA-Wilson-Prostate-Primaries-Not-Spread_n=21_Norm-MiSeq.xls")[8,c(-1,-23)])
nospread_p_ids_sheet=sapply(strsplit(nospread_p_ids_sheet,"[-_]"),function(x)paste0(x[2],"-",x[3]))


###MIRNA NAMES
##2102 miRNAs probes all together between the two, they match and its ok

mirnas.names.1<-lapply(d,function(x){
  temp<-read_xls(x,skip=8,sheet=2)
  unlist(temp[-1:-2,1])
  
})


mirnas.names.2<-lapply(Sys.glob("Data/data_2/*"),function(x){
  temp<-read_xls(x,skip=8,sheet=2)
  unlist(temp[-1:-2,1])
  
})

#sum(mirnas.names.1[[1]]==mirnas.names.2[[1]])







####Load in all the IDs and names, etc for each sample
ids.1<-read_xlsx("Data/sample_info/sample_info.xlsx") %>% 
  mutate(ID1=sapply(strsplit(S.number,"[.]"),function(x)paste(x[1:2],collapse="."))) %>% 
  mutate(location=gsub("[0-9]","",location)) %>% 
  mutate(tissue=ifelse(tissue=="lymph_node","Lymph Node",tissue)) %>%
  mutate(tissue=ifelse(tissue=="primary","Primary",ifelse(tissue=="secondary","Secondary",tissue))) %>%
  mutate(location=ifelse(location=="primary","Primary",ifelse(location=="secondary","Secondary",location))) %>%
  mutate(tissue1=ifelse(tissue=="Primary"|tissue=="Secondary","Tumor","Lymph Node")) %>% 
  mutate(location1=ifelse(location2=="primary"|location2=="Secondary","Tumor",location2)) %>%
  select(ID1,ID2=S.number,tissue,tissue1,location,location1)
#grepl("[0-9]",str_sub(location,-1


ids.2<-read_xlsx("Data/sample_info/B2_SampleList.xlsx") %>% na.omit()  
names(ids.2)<-c("Sample_ID","Block_ID","SA","MRN","Slides_num","Location")
ids.2<-ids.2 %>%
  mutate(temptumor=sapply(strsplit(`Block_ID`,"[-]"),function(x)x[length(x)])) %>%
  mutate(temptumor=gsub(" ","",temptumor)) %>%
  mutate(temptumor=toupper(temptumor)) %>%
  mutate(temptumor=ifelse(Location=='Secondary',"Secondary",temptumor)) %>%
  mutate(temptumor=ifelse(temptumor=='RP',"Primary",temptumor)) %>%
  mutate(temptumor=ifelse(temptumor %in% c("B1","D1","G1"),"Normal",temptumor))  %>%
  separate('Sample_ID',into=c("TEMP1","TEMP2"),sep="-") %>% 
  separate('Block_ID',into=c("B1","B2"),sep="[- ]") %>% 
  mutate(ID2=paste0(TEMP1,".",TEMP2,".",B1))  %>%
  mutate(ID1=paste0(TEMP1,".",TEMP2))  %>%
  mutate(tissue=ifelse(temptumor %in% c("Primary","Secondary"),temptumor,"Lymph Node")) %>%
  mutate(tissue1=ifelse(temptumor %in% c("Primary","Secondary"),"Tumor","Lymph Node")) %>%
  mutate(location1=ifelse(temptumor %in% c("Primary","Secondary"),"Tumor",temptumor)) %>%
  select(ID1,ID2,tissue,tissue1,location=temptumor,location1)


###Ids for Normal
normal_p_ids=read_xlsx("Data/data_3/Simple Prostatectomy Patients Info. 1029.xlsx")




##Read in covariate data
covar_plnd<-read_xlsx("Data/sample_info/PLND for Darryl.xlsx") %>% 
  filter(!is.na(`Accession ID`)) %>% 
  select(ID=`Accession ID`,Age=`Age at Surgery`,PSA=PrePSA,PSA_cat='Post-op PSA category',TimeBCR='Time to BCR (mos)',Tumorsize=`Dominant Tumor Size (cm)`,LN_num="#LN",LN_pro=`Lymph node`) %>%
  mutate(LN_pro=sapply(strsplit(LN_pro,"/"),function(x)as.numeric(x[1])/as.numeric(x[2])))


###BIND ALL B1
mirnas.b1<-do.call(cbind,mirnas_b1)
mirnas.b1<-mirnas.b1[,names(mirnas.b1) %in% ids.1$ID2]
mirnas.b1<-mirnas.b1 %>% select(ids.1$ID2)



names(mirnas.b2)<-ids.2$ID2




##Final ID
ids=bind_rows(ids.1,ids.2)



ids=ids %>% group_by(ID1) %>% mutate(ID3=ifelse(any(tissue %in% "Secondary"),"MF","SF")) %>% ungroup() %>% group_by(ID3) %>%
  mutate(T1=as.numeric(factor(ID1,levels=unique(ID1)))) %>% ungroup() %>% mutate(ID3=paste0(ID3,T1)) %>% select(-T1) %>%
  mutate(ID4=make.unique(paste0(ID3," ",location),sep=" "))




normal_ids=tibble(ID1=normal_p_ids_sheet,ID2=normal_p_ids_sheet,tissue="Normal",tissue1="Normal",
                  location="Normal",location1="Normal",
                  ID3=normal_p_ids_sheet,ID4=normal_p_ids_sheet)

nospread_ids=tibble(ID1=nospread_p_ids_sheet,ID2=nospread_p_ids_sheet,tissue="PCa_Nospread",tissue1="PCa_Nospread",
                    location="PCa_Nospread",location1="PCa_Nospread",
                    ID3=nospread_p_ids_sheet,ID4=nospread_p_ids_sheet)


ids=bind_rows(ids,normal_ids,nospread_ids)




###Filter out the  Control Probes and low Performing
###normalize first and then extract the samples that we are interested in 
##11 probes are removed for all leaving 2083
##updated analysis in 2020-only 3 are removed due 

mirnas.dt<-cbind(mirnas.b1,mirnas.b2,normal_p,nospread_p)
rownames(mirnas.dt)=mirnas.names.1[[1]]
colnames(mirnas.dt)=ids$ID2


#write_rds(mirnas.dt,"Data/Rdata/mirnas_ALL_dt.rds",compress="xz")



mirnas.dt=mirnas.dt[-1:-19,]

rmvec<-!apply(mirnas.dt,1,sum)<10

mirnas.dt1<-mirnas.dt[rmvec,] 



##Also remove names
mirnas.names<-mirnas.names.1[[1]][-1:-19]
mirnas.names<-mirnas.names[rmvec]


##Add rownames
rownames(mirnas.dt1)=mirnas.names

##Add col names
names(mirnas.dt1)=ids$ID2



##Scale all individuals to CPM and take the log2 base of all the samples
sum_tot=apply(mirnas.dt1,2,sum)
mirnas.dt2<-apply(mirnas.dt1,2,function(x){
  x1=((x*1e6)/sum_tot)
  x1=ifelse(x1<2,2,x1)
  log(x1,base=2)
})



###Run GGPLOT for all of the  batches to see 
#sapply(mirnas_b1,function(x)sum(names(x) %in% ids$ID2))

batch_counts=data.frame(Counts=apply(mirnas.dt1,2,sum),
                        Batch=c(rep("B1",14),
                                rep("B2",8),rep('B3',11),
                                rep(c('B4','B5','B6'),each=23),
                                rep("B7",24),
                                rep("B8",23),
                                rep("B9","21")))


#ggplot(batch_counts,aes(x=Batch,y=Counts))+
#  geom_violin() +
#  geom_point(position=position_jitter(width=0.2),size=0.3) +
#  ylab("Total Counts")




covar<-data.frame(ids,batch_counts) %>% 
  mutate(ID5=gsub("[.]","-",ID1)) %>% 
  left_join(.,covar_plnd,by=c("ID5"="ID"))






###For Plotting purposes 
mirnas.de<-mirnas.dt1

mirnas_plot<-DESeqDataSetFromMatrix(countData=mirnas.de,
                            colData=covar,
                            design=~Batch)


de.vst<-assay(vst(mirnas_plot,fitType = "local"))
de.vst.raw=vst(mirnas_plot,fitType = "local")
DESeq2::plotPCA(de.vst.raw,intgroup="Batch",ntop=1000,returnData=T)


#head(as.data.frame(x))

###CHECKING BATCH EFFECTS
modcombat = model.matrix(~1,data = covar)
combat_mirna = ComBat(dat=de.vst, batch=covar$Batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)


mirnas_de_pc=prcomp(t(mirnas.de))
combat_pc=prcomp(t(combat_mirna))

ggplot(data.frame(mirnas_de_pc$x)) + geom_point(aes(x=PC1,y=PC2)) 
ggplot(data.frame(combat_pc$x)) + geom_point(aes(x=PC1,y=PC2)) 



edata_com100<-combat_mirna[order(-apply(combat_mirna,1,sd)),]
edata_com100<-edata_com100[1:100,]


write_rds(ids,"Data/Rdata/ids.rds",compress = "xz")
write_rds(batch_counts,"Data/Rdata/batch_counts.Rds",compress = "xz")
write_rds(covar,"Data/Rdata/covar.Rds",compress = "xz")
write_rds(mirnas.de,"Data/Rdata/mirnas.de.Rds","xz")
write_rds(combat_mirna,"Data/Rdata/combat_mirna.Rds","xz")

