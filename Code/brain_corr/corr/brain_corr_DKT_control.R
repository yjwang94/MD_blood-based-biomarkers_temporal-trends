#提取脑区数据
library(stringr)
fpath = "/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/Neuroimaging/"
volume = read.csv(paste(fpath,"Freesurfer_DKT_bl_data.csv",sep = ""))
area1 = read.csv(paste(fpath,"mddadult_case-controls_CortThick.csv",sep = ""))
DKT_field = read.csv(paste(fpath,"Freesurfer DKT.csv",sep = ""))
DKT_field_decription = read.csv(paste(fpath,"Freesurfer DKT1.csv",sep = ""))
DKT_field_volume = DKT_field[c(125:186,230:283),]
area1_lst = area1$Structure
for (i in c(1:68)) {
  area1_lst[i] = substr(area1[i,2],3,nchar(area1[i,2]))
}

DKT_field_lst = DKT_field$Description[c(125:186,230:283)]
for (i in c(1:116)) {
  DKT_field_lst[i] = str_split(DKT_field_lst[i],' ',simplify = T)[1,3]
}

DKT_field_volume1 = DKT_field_volume[DKT_field_lst %in% area1_lst,]
DKT_field_volume2 = DKT_field[c(233:236,247:248,258:259,262:263,265:268,272:273),]
DKT_field_volume_area1 = DKT_field_volume$Field.ID[DKT_field_lst %in% area1_lst]
DKT_field_volume_area2 = DKT_field$Field.ID[c(233:236,247:248,258:259,262:263,265:268,272:273)]

volume_field_lst = colnames(volume)[2:286]
for (i in c(1:285)) {
  volume_field_lst[i] = substr(volume_field_lst[i],2,6)
}
TIV = volume[colnames(volume) == "x26521_2_0"]
colnames(TIV) = "TIV"

volume_DKT = volume[,2:286][,volume_field_lst %in% DKT_field_volume_area1 | volume_field_lst %in% DKT_field_volume_area2]


volume_DKT_lst = colnames(volume_DKT)
for (i in 1:78) {
  volume_DKT_lst[i] = substr(volume_DKT_lst[i],2,6)
}
colnames(volume_DKT) = volume_DKT_lst


volume_DKT_lst = colnames(volume_DKT)
for (i in 1:78) {
  volume_DKT_lst[i] = DKT_field_decription$Description1[DKT_field_decription$Field.ID == volume_DKT_lst[i]]
}
colnames(volume_DKT) = volume_DKT_lst
volume_DKT = cbind(volume[,1], volume_DKT)
colnames(volume_DKT)[1] = "eid"

volume_DKT = cbind(volume_DKT, TIV)

#计算左右脑平均
#substr(colnames(volume_DKT)[2],-1,-2)
#str_split(colnames(volume_DKT)[2],'_',simplify = T)[1,2]
volume_DKT_mean = volume_DKT[1]
volume_DKT_mean1 = data.frame(matrix(nrow = 43173,ncol = 8))
for (i in 1:8) {
  volume_DKT_mean1[i]=(volume_DKT[i+1]+volume_DKT[i+9])/2
}
colnames(volume_DKT_mean1)=c("Lateral-Ventricle", "Thalamus-Proper", "Caudate",                 
                       "Putamen", "Pallidum","Hippocampus", "Amygdala",                
                       "Accumbens-area")
volume_DKT_mean2 = data.frame(matrix(nrow = 43173,ncol = 31))
for (i in 1:31) {
  volume_DKT_mean2[i]=(volume_DKT[i+17]+volume_DKT[i+48])/2
}
colnames(volume_DKT_mean2)=c("caudalan terior cingulate","caudal middle frontal",
                             "cuneus", "entorhinal", "fusiform",                
                             "inferior parietal", "inferior temporal",        
                             "isthmus cingulate", "lateral occipital",        
                             "lateral orbitofrontal", "lingual",                 
                             "medial orbitofrontal", "middle temporal",          
                             "para hippocampal", "para central",             
                             "pars opercularis", "pars orbitalis",           
                             "pars triangularis", "perical carine",           
                             "post central", "posterior cingulate",      
                             "pre central", "pre cuneus",               
                             "rostral anterior cingulate", "rostral middle frontal",    
                             "superior frontal", "superior parietal",      
                             "superior temporal", "supramarginal",           
                             "transverse temporal", "insula")

volume_DKT_mean = cbind(volume_DKT_mean,volume_DKT_mean1,volume_DKT_mean2,TIV)


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

#MD显著变量的并集
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
control_sig_volume_DKT_mean_cov = merge(control_blood, volume_DKT_mean, on = "eid")

#回归协变量
#eliminate the influence of covariance
data_cov_pheno = control_sig_volume_DKT_mean_cov
data_cov_pheno_elimate = control_sig_volume_DKT_mean_cov
for (i in (47:(46+length(all_significant_variable_sorted_lst)))) {
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
for (i in ((47+length(all_significant_variable_sorted_lst)):(length(data_cov_pheno)-1))){
  glm1 = glm(data_cov_pheno[,i] ~ Gender + TD_Index + Age + SMK_Status1 + SMK_Status2 + BMI +
               ALC_Status1 + ALC_Status2 + Education3 + Education4 + 
               Education5 + Education7 + Education8 + TIV +
               Ethnicity1 + Ethnicity2 + Ethnicity3 + Sites11001 + Sites11002 + Sites11003 +
               Sites11004 + Sites11005 + Sites11006 + Sites11007 + Sites11008 + Sites11009 + 
               Sites11010 + Sites11011 + Sites11012 + Sites11013 + Sites11014 + Sites11016 + 
               Sites11017 + Sites11018 + Sites11020 + Sites11021 + Sites11022 + Sites11023, data = data_cov_pheno)
  
  res = residuals(glm1)
  
  new_data = res+glm1$coefficients[1]
  
  data_cov_pheno_elimate[complete.cases(data_cov_pheno_elimate[,i]),][,i] = new_data
}
#calculate zscores
zdata = data_cov_pheno_elimate
for (i in (47:(46+length(all_significant_variable_sorted_lst)))) {
  zdata[complete.cases(zdata[,i]),][,i] = (zdata[complete.cases(zdata[,i]),][,i] - mean(zdata[complete.cases(zdata[,i]),][,i]))/sd(zdata[complete.cases(zdata[,i]),][,i])
}
#colnames(zdata)[2] = 'duration'
MD_pheno_sig_volume_DKT_mean_cov_elimate = zdata
#计算zscore

#画热图
library(reshape2)
library(dplyr)
library(ggplot2)

corrs_data = data.frame(matrix(nrow = 39))
pvalues_data = data.frame(matrix(nrow = 39))
for (i in c(47:(46+length(all_significant_variable_sorted_lst)))) {
  corrs = data.frame()
  pvalues = data.frame()
  for (j in c((47+length(all_significant_variable_sorted_lst)):(length(all_significant_variable_sorted_lst)+85))) {
    data_complete = MD_pheno_sig_volume_DKT_mean_cov_elimate[complete.cases(MD_pheno_sig_volume_DKT_mean_cov_elimate[i]),]
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
colnames(corrs_data) = colnames(MD_pheno_sig_volume_DKT_mean_cov_elimate)[47:(46+length(all_significant_variable_sorted_lst))]
rownames(corrs_data) = colnames(MD_pheno_sig_volume_DKT_mean_cov_elimate)[(47+length(all_significant_variable_sorted_lst)):(length(all_significant_variable_sorted_lst)+85)]
colnames(pvalues_data) = colnames(MD_pheno_sig_volume_DKT_mean_cov_elimate)[47:(46+length(all_significant_variable_sorted_lst))]
rownames(pvalues_data) = colnames(MD_pheno_sig_volume_DKT_mean_cov_elimate)[(47+length(all_significant_variable_sorted_lst)):(length(all_significant_variable_sorted_lst)+85)]
rownames(corrs_data)[2] = "Thalamus"
rownames(corrs_data)[8] = "Accumbens"
library(Hmisc)
rownames(corrs_data) = capitalize(rownames(corrs_data))
rownames(pvalues_data) = capitalize(rownames(pvalues_data))

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
ggsave("DKT_control.pdf",width = 20,height = 10)


