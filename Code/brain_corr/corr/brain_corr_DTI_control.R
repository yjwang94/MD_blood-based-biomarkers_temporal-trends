#DTI
library(stringr)
library(dplyr)
fpath = "/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/Neuroimaging/"
DTI_field = read.csv(paste(fpath,"Diffusion_MRI_weighted_means.csv",sep = ""))
DTI = read.csv(paste(fpath,"Diffusion_MRI_bl_data.csv",sep = ""))
tmp = colnames(DTI)
for (i in c(2:244)) {
  tmp[i] = sub('x', '', tmp[i])
  tmp[i] = substr(tmp[i], 1, 5)
}
colnames(DTI) = tmp

DTI_field_FA = DTI_field[c(1:27),]
DTI_field_MD = DTI_field[c(163:189),]
DTI_field_ICVF = DTI_field[c(28:54),]

#添加缩写
for (i in c(1:27)) {
  a = substr(DTI_field_FA$Description[i],27,nchar(DTI_field_FA$Description[i]))
  b = str_split(a,' ',simplify = T)
  b<-unlist(b)
  c = list(substr(b, start = 1,stop = 1))
  d = unlist(c)
  e = ""
  for (j in c(1:length(d)-1)) {
    e = paste(e,d[j],sep = "")
  }
  e = toupper(e)
  e = paste("FA",e)
  DTI_field_FA[i,3] = e
}
colnames(DTI_field_FA)[3] = c("abbreviation")
DTI_field_FA[9,3] = "FA FMA"
DTI_field_FA[10,3] = "FA FMI"

for (i in c(1:27)) {
  a = substr(DTI_field_MD$Description[i],27,nchar(DTI_field_MD$Description[i]))
  b = str_split(a,' ',simplify = T)
  b<-unlist(b)
  c = list(substr(b, start = 1,stop = 1))
  d = unlist(c)
  e = ""
  for (j in c(1:length(d)-1)) {
    e = paste(e,d[j],sep = "")
  }
  e = toupper(e)
  e = paste("MD",e)
  DTI_field_MD[i,3] = e
}
colnames(DTI_field_MD)[3] = c("abbreviation")
DTI_field_MD[9,3] = "MD FMA"
DTI_field_MD[10,3] = "MD FMI"


for (i in c(1:27)) {
  a = substr(DTI_field_ICVF$Description[i],29,nchar(DTI_field_ICVF$Description[i]))
  b = str_split(a,' ',simplify = T)
  b<-unlist(b)
  c = list(substr(b, start = 1,stop = 1))
  d = unlist(c)
  e = ""
  for (j in c(1:length(d)-1)) {
    e = paste(e,d[j],sep = "")
  }
  e = toupper(e)
  e = paste("ICVF",e)
  DTI_field_ICVF[i,3] = e
}
colnames(DTI_field_ICVF)[3] = c("abbreviation")
DTI_field_ICVF[9,3] = "ICVF FMA"
DTI_field_ICVF[10,3] = "ICVF FMI"
DTI_field_all = rbind(DTI_field_FA, DTI_field_MD, DTI_field_ICVF)

#计算左右脑平均
unique3 = c("25498","25499","25504","25525","25526","25531","25660","25661","25666")
DTI_field_all1 = DTI_field_all[!DTI_field_all$Field.ID %in% unique3,]
DTI_all = DTI[,colnames(DTI) %in% DTI_field_all$Field.ID & !(colnames(DTI) %in% unique3)]
DTI_3 = DTI[,colnames(DTI) %in% unique3]
colnames(DTI_3) = c("FA FMA","FA FMI","FA MC","MD FMA","MD FMI", "MD MC", "ICVF FMA","ICVF FMI","ICVF MC")
DTI_mean1 = data.frame(matrix(nrow = 40747,ncol = 39))

for (i in seq(1,72,2)) {
  DTI_mean1[(i+1)/2]=(DTI_all[i]+DTI_all[i+1])/2
}
colnames(DTI_mean1) = unique(DTI_field_all1$abbreviation)
DTI_mean1 = cbind(DTI_mean1,DTI_3)
DTI_mean1 = DTI_mean1[, order(colnames(DTI_mean1))]
DTI_mean = cbind(DTI[1],DTI_mean1)
DTI_mean = DTI_mean[c(1:46)]

#匹配Serum & NMR
UKB_pheno = read.csv("/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/PhenoData2Analysis/UKB_pheno.csv")
UKB_blood_dict = read.csv("/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/PhenoData2Analysis/UKB_Blood_dict.csv")
NMR_dict = read.csv("/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/PhenoData2Analysis/abbreviation_NMR.csv")
Serum_dict = read.csv("/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/PhenoData2Analysis/Serum_dict1.csv")
UKB_blood_dict$FieldID = paste(UKB_blood_dict$FieldID, '-0.0',sep = '')
Serum_dict$FieldID = paste(Serum_dict$FieldID,'-0.0', sep = '')
UKB_blood_dict$Field = tolower(UKB_blood_dict$Field)
NMR_dict$Description = tolower(NMR_dict$Description)
Serum_dict$Field = tolower(Serum_dict$Field)

# rename pheno of UKB_pheno
tmp = colnames(UKB_pheno)[2:383]
for (i in (1:382)) {
  tmp[i] = sub('X', '', tmp[i])
  tmp[i] = sub('\\.', '\\-', tmp[i])
}
colnames(UKB_pheno)[2:383] = tmp

#将NMR_dict加上FieldID
for (i in NMR_dict$Description) {
  j = UKB_blood_dict$FieldID[UKB_blood_dict$Field == i]
  NMR_dict$FieldID[NMR_dict$Description == i] = j
}

data_NMR = UKB_blood_dict[UKB_blood_dict$Main_Category == "Metabolomics",]
for (i in data_NMR$Field) {
  j = NMR_dict$Abbreviation[NMR_dict$Description == i]
  UKB_blood_dict$field[UKB_blood_dict$Field == i] = j
}


#Serum_dict$Sub_Category1 = tolower(Serum_dict$Sub_Category1)

data_Serum = UKB_blood_dict[UKB_blood_dict$Main_Category == "Serum",]
for (i in data_Serum$FieldID) {
  j = Serum_dict$Sub_Category1[Serum_dict$FieldID == i]
  UKB_blood_dict$Sub_Category[UKB_blood_dict$FieldID == i] = j
}
for (i in data_Serum$FieldID) {
  j = Serum_dict$Abbreviation[Serum_dict$FieldID == i]
  UKB_blood_dict$field[UKB_blood_dict$FieldID == i] = j
}


Serum_Field = Serum_dict[,c(1, 6)]
NMR_Field = NMR_dict[,c(8,1)]
data_FieldID = rbind(Serum_Field, NMR_Field)
UKB_pheno227 = UKB_pheno[, c('eid', UKB_blood_dict$FieldID)]

fpath = "/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/"
S0_TargetData = read.csv(paste(fpath,"raw_data/S0_TargetData.csv",sep = ""))
cov = read.csv(paste(fpath,"raw_data/Covariates_Data_Imputed.csv",sep = ""))
data_cov = cov
data_cov$treat = 0
MD_list = c("AX","BIP","DEP","SCH")
for (MD in MD_list) {
  S0_TargetData1 = S0_TargetData[S0_TargetData[paste(MD,"_date",sep="")] != "",]
  data_cov$treat[data_cov$eid %in% S0_TargetData1$eid] = 1
}
#对协变量进行处理
data_cov$SMK_Status = as.character(data_cov$SMK_Status)
data_cov$ALC_Status = as.character(data_cov$ALC_Status)
data_cov$Gender = as.character(data_cov$Gender)
data_cov$Education = as.character(data_cov$Education)
data_cov$Ethnicity = as.character(data_cov$Ethnicity)
data_cov$Sites = as.character(data_cov$Sites)

dumy <- model.matrix( ~ SMK_Status - 1, data = data_cov)
data_cov <- cbind(data_cov,dumy[,-1])
dumy <- model.matrix( ~ ALC_Status - 1, data = data_cov)
data_cov <- cbind(data_cov,dumy[,-1])
dumy <- model.matrix( ~ Education - 1, data = data_cov)
data_cov <- cbind(data_cov,dumy[,-1])
dumy <- model.matrix( ~ Ethnicity - 1, data = data_cov)
data_cov <- cbind(data_cov,dumy[,-1])
dumy <- model.matrix( ~ Sites - 1, data = data_cov)
data_cov <- cbind(data_cov,dumy[,-1])
#连续变量规范化处理
continues = c("TD_Index", "Age", "BMI")

#data_cov1 = data_cov
for (i in continues) {
  data_cov[,i] = (data_cov[,i]-mean(data_cov[,i])) / sd(data_cov[,i])
}

case = data_cov[data_cov$treat == 1,]
control = data_cov[data_cov$treat == 0,] 
Blood_immune = read.csv(paste(fpath,'raw_data/disease_need_exclusion_data/','Blood and immune disorders.csv',sep = ""))
Blood_immune_eid_lst = Blood_immune$eid[Blood_immune$target_y == 1 & Blood_immune$BL2Target_yrs <0]
Mental_behavioural = read.csv(paste(fpath,'raw_data/disease_need_exclusion_data/','Mental and behavioural disorders.csv',sep = ""))
Mental_behavioural_eid_lst = Mental_behavioural$eid[Mental_behavioural$target_y == 1 & Mental_behavioural$BL2Target_yrs <0]
Infections = read.csv(paste(fpath,'raw_data/disease_need_exclusion_data/','Diseases of Infections.csv',sep = ""))
Infections_eid_lst = Infections$eid[Infections$target_y == 1 & Infections$BL2Target_yrs <0]
Nervous = read.csv(paste(fpath,'raw_data/disease_need_exclusion_data/','Nervous system disorders.csv',sep = ""))
Nervous_eid_lst = Nervous$eid[Nervous$target_y == 1 & Nervous$BL2Target_yrs < 0]
other_disease = unique(c(Blood_immune_eid_lst,Mental_behavioural_eid_lst,Infections_eid_lst,Nervous_eid_lst ))
control_noohter = control[!control$eid %in% other_disease,]

all_significant_variable_sorted_lst = c()
for (MD in MD_list) {
  significant_variable = read.csv(paste("/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/CaseControl/S6_cluster_clustra/S62_cluster_data_",MD,".csv",sep = ""))
  significant_variable_sorted <- arrange(significant_variable, cluster)
  significant_variable_lst = significant_variable$X
  significant_variable_sorted_lst = significant_variable_sorted$X
  all_significant_variable_sorted_lst = unique(append(all_significant_variable_sorted_lst, significant_variable_sorted_lst))
}
control_blood = merge(control_noohter, UKB_pheno227[,c("eid",all_significant_variable_sorted_lst)], by = 'eid')
#把脑区体积和control变量数据进行匹配按照eid
control_sig_volume_DTI_mean_cov = merge(control_blood, DTI_mean, on = "eid")
#回归协变量
#eliminate the influence of covariance
data_cov_pheno = control_sig_volume_DTI_mean_cov
data_cov_pheno_elimate = control_sig_volume_DTI_mean_cov
for (i in (47:length(data_cov_pheno))){
  glm1 = glm(data_cov_pheno[,i] ~ Gender + TD_Index + Age + SMK_Status1 + SMK_Status2 + BMI +
               ALC_Status1 + ALC_Status2 + Education3 + Education4 + 
               Education5 + Education7 + Education8 + 
               Ethnicity1 + Ethnicity2 + Ethnicity3 + Sites11001 + Sites11002 + Sites11003 +
               Sites11004 + Sites11005 + Sites11006 + Sites11007 + Sites11008 + Sites11009 + 
               Sites11010 + Sites11011 + Sites11012 + Sites11013 + Sites11014 + Sites11016 + 
               Sites11017 + Sites11018 + Sites11020 + Sites11021 + Sites11022 + Sites11023, data = data_cov_pheno)
  
  res = residuals(glm1)
  
  new_data = res+glm1$coefficients[1]
  
  data_cov_pheno_elimate[complete.cases(data_cov_pheno_elimate[,i]),][,i] = new_data
} #take so long time

#calculate zscores
zdata = data_cov_pheno_elimate
for (i in (47:(46+length(all_significant_variable_sorted_lst)))) {
  zdata[complete.cases(zdata[,i]),][,i] = (zdata[complete.cases(zdata[,i]),][,i] - mean(zdata[complete.cases(zdata[,i]),][,i]))/sd(zdata[complete.cases(zdata[,i]),][,i])
}
#colnames(zdata)[2] = 'duration'
MD_pheno_sig_volume_DTI_mean_cov_elimate = zdata
#计算zscore

#画热图
library(reshape2)
library(dplyr)
library(ggplot2)

corrs_data = data.frame(matrix(nrow = 45))
pvalues_data = data.frame(matrix(nrow = 45))
for (i in c(47:(46+length(all_significant_variable_sorted_lst)))) {
  corrs = data.frame()
  pvalues = data.frame()
  for (j in c((47+length(all_significant_variable_sorted_lst)):(length(all_significant_variable_sorted_lst)+91))) {
    data_complete = MD_pheno_sig_volume_DTI_mean_cov_elimate[complete.cases(MD_pheno_sig_volume_DTI_mean_cov_elimate[i]),]
    corr = cor.test(data_complete[,i],data_complete[,j],method = "spearman")$estimate
    pvalue = cor.test(data_complete[,i],data_complete[,j],method = "spearman")$p.value
    corrs = rbind(corrs,corr)
    pvalues = rbind(pvalues,pvalue)
  }
  corrs_data = cbind(corrs_data, corrs)
  pvalues_data = cbind(pvalues_data, pvalues)
}
corrs_data = corrs_data[,-1]
pvalues_data = pvalues_data[,-1]
colnames(corrs_data) = colnames(MD_pheno_sig_volume_DTI_mean_cov_elimate)[47:(46+length(all_significant_variable_sorted_lst))]
rownames(corrs_data) = colnames(MD_pheno_sig_volume_DTI_mean_cov_elimate)[(47+length(all_significant_variable_sorted_lst)):(length(all_significant_variable_sorted_lst)+91)]
colnames(pvalues_data) = colnames(MD_pheno_sig_volume_DTI_mean_cov_elimate)[47:(46+length(all_significant_variable_sorted_lst))]
rownames(pvalues_data) = colnames(MD_pheno_sig_volume_DTI_mean_cov_elimate)[(47+length(all_significant_variable_sorted_lst)):(length(all_significant_variable_sorted_lst)+91)]
#rownames(corrs_data)[2] = "Thalamus"
#rownames(corrs_data)[8] = "Accumbens"
library(Hmisc)
#rownames(corrs_data) = capitalize(rownames(corrs_data))
#rownames(pvalues_data) = capitalize(rownames(pvalues_data))

#对FieldID进行替换成description Serum NMR
tmp = colnames(corrs_data)
for (i in c(1:58)) {
  if(tmp[i] %in% NMR_Field$FieldID){
    tmp[i] = NMR_Field$Abbreviation[NMR_Field$FieldID == tmp[i]]
  }
  if(tmp[i] %in% Serum_Field$FieldID){
    tmp[i] = Serum_Field$Abbreviation[Serum_Field$FieldID == tmp[i]]
  }
}
colnames(corrs_data) = tmp
colnames(pvalues_data) = tmp

pvalues_data_bfi = pvalues_data
for (i in c(1:39)) {
  pvalues_data_bfi[i,] = p.adjust(pvalues_data[i,],
                                  method = "bonferroni")
}

pvalues_matrix = as.matrix(pvalues_data_bfi)
corrs_matrix = as.matrix(corrs_data)

p <- melt(pvalues_matrix)
d <- melt(corrs_matrix)

d <- cbind(d,p[3])
colnames(d) = c("Var1","Var2","value","value1")
d <- d %>%
  mutate(text = case_when(
    0.005 < value1 & value1 < 0.05 ~ paste("*"), 
    0.0005 < value1 & value1 <= 0.005 ~ paste("**"),
    value1 <= 0.0005 ~ paste("***")))
d <- d[,-4]
class(d)
ggplot(d, aes(Var2, Var1)) + 
  geom_tile(aes(fill = value), colour = "grey", size = 0.5)+
  scale_fill_gradient2(low = "#0f86a9",mid = "white",high = "#FC8452") +
  geom_text(aes(label=text),col ="black",size = 5,vjust = 0.8) +
  theme_minimal() + 
  theme(axis.title.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.title.y=element_blank(), 
        axis.text.x = element_text(hjust = 1, size = 15, face = "bold",angle=45), 
        axis.text.y = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold")) + 
  labs(fill =paste0(" * p_bfi < 0.05","\n\n","** p_bfi < 0.005","\n\n",
                    "*** p_bfi < 0.0005","\n\n","corr"))    

setwd("/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/Figures/brain_blood_corr/corr/")
ggsave(paste("DTI_control.pdf",sep = ""),width = 15,height = 10)
