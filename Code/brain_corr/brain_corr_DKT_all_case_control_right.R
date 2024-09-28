#提取脑区数据
library(stringr)
library(reshape2)
library(dplyr)
library(ggplot2)
#install.packages("htmltools")
library(Hmisc)

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

#############
UKB_pheno = read.csv("/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/PhenoData2Analysis/UKB_pheno.csv")
UKB_blood_dict = read.csv("/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/PhenoData2Analysis/UKB_Blood_dict.csv")
NMR_dict = read.csv("/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/PhenoData2Analysis/abbreviation_NMR.csv")
Serum_dict = read.csv("/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/PhenoData2Analysis/Serum_dict1.csv")
#UKB_blood_dict$FieldID = paste(UKB_blood_dict$FieldID, '-0.0',sep = '')
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

#MD血液指标和脑区体积进行配对， with 39个脑区体积 +协变量
#并且剔除其他疾病的人群53903
MD_list = c("AX","DEP","BIP","SCH")
MD = "BIP"
betas_all = data.frame(matrix(nrow = 39))
pvalues_all = data.frame(matrix(nrow = 39))
low_CI_all = data.frame(matrix(nrow = 39))
up_CI_all = data.frame(matrix(nrow = 39))
for (MD in MD_list) {
  #原始数据而非zscore
  MD_zscore = read.csv(paste("/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/CaseControl/",MD,"/S2_data_cov_pheno_forPvalue.csv",sep = ""))
  disease_exclude_eid = read.csv("/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/Neuroimaging/disease exclude eid.csv")
  disease_exclude_eid_lst = disease_exclude_eid$eid
  MD_zscore = MD_zscore[!(MD_zscore$eid %in% disease_exclude_eid_lst),]
  #协变量处理
  data_cov = MD_zscore
  
  #dummy 分类变量
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
  data_cov1 = cbind(data_cov[,1:14],data_cov[,242:274],data_cov[,15:241])
  tmp = colnames(data_cov1)[48:274]
  for (i in (1:227)) {
    tmp[i] = sub('X', '', tmp[i])
    tmp[i] = sub('\\.', '\\-', tmp[i])
  }
  colnames(data_cov1)[48:274] = tmp
  
  MD_zscore = data_cov1
  
  significant_variable = read.csv(paste("/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/CaseControl/S6_cluster_clustra/S62_cluster_data_",MD,".csv",sep = ""))
  significant_variable_sorted <- arrange(significant_variable, cluster)
  significant_variable_lst = significant_variable$X
  significant_variable_sorted_lst = significant_variable_sorted$X
  # rename pheno of UKB_pheno
  colnames(MD_zscore)[48:(47+length(significant_variable_lst))] = significant_variable_lst
  #对MD_zscore中的blood变量按照cluster排序
  MD_zscore_arrange = cbind(MD_zscore[,1:47],MD_zscore[,significant_variable_sorted_lst])
  #把脑区体积和变量数据进行匹配按照eid
  MD_pheno_sig_volume_DKT_mean_cov = merge(MD_zscore_arrange, volume_DKT_mean, on = "eid")
  
  
  #拟合线性模型
  betas_data = data.frame(matrix(nrow = 39))
  pvalues_data = data.frame(matrix(nrow = 39))
  low_CI_data = data.frame(matrix(nrow = 39))
  up_CI_data = data.frame(matrix(nrow = 39))
  i=48
  j=79
  for (i in (48:(47+length(significant_variable_lst)))){#blood biomarkers
    betas = data.frame()
    pvalues = data.frame()
    low_CIs = data.frame()
    up_CIs = data.frame()
    for (j in ((48+length(significant_variable_lst)):(length(significant_variable_lst)+86))) {#脑区
      complete = MD_pheno_sig_volume_DKT_mean_cov[complete.cases(MD_pheno_sig_volume_DKT_mean_cov[,i]),]
      lm1 = lm(complete[,j] ~ complete[,'group']*complete[,i]+ Gender + TD_Index + Age + SMK_Status1 + SMK_Status2 + TIV +
                 ALC_Status1 + ALC_Status2 + Education3 + Education4 + 
                 Education5 + Education7 + Education8 + 
                 Ethnicity1 + Ethnicity2 + Ethnicity3 + Sites11001 + Sites11002 + Sites11003 +
                 Sites11004 + Sites11005 + Sites11006 + Sites11007 + Sites11008 + Sites11009 +
                 Sites11010 + Sites11011 + Sites11012 + Sites11013 + Sites11014 + Sites11016 + 
                 Sites11017 + Sites11018 + Sites11020 + Sites11021 + Sites11022 + Sites11023, data = complete)
      beta = lm1$coefficients[length(lm1$coefficients)]
      sumres <- summary(lm1)
      pvalue = sumres$coefficients[length(sumres$coefficients)]
      low_CI = confint(lm1)[length(confint(lm1))/2]
      up_CI = confint(lm1)[length(confint(lm1))]
      betas = rbind(betas,beta)
      pvalues = rbind(pvalues,pvalue)
      low_CIs = rbind(low_CIs, low_CI)
      up_CIs = rbind(up_CIs, up_CI)
    }
    betas_data = cbind(betas_data, betas)
    pvalues_data = cbind(pvalues_data, pvalues)
    low_CI_data = cbind(low_CI_data, low_CIs)
    up_CI_data = cbind(up_CI_data, up_CIs)
  }
  
  betas_data = betas_data[,-1]
  pvalues_data = pvalues_data[,-1]
  low_CI_data = low_CI_data[,-1]
  up_CI_data = up_CI_data[,-1]
  colnames(betas_data) = colnames(MD_pheno_sig_volume_DKT_mean_cov)[48:(47+length(significant_variable_lst))]
  rownames(betas_data) = colnames(MD_pheno_sig_volume_DKT_mean_cov)[(48+length(significant_variable_lst)):(length(significant_variable_lst)+86)]
  colnames(pvalues_data) = colnames(MD_pheno_sig_volume_DKT_mean_cov)[48:(47+length(significant_variable_lst))]
  rownames(pvalues_data) = colnames(MD_pheno_sig_volume_DKT_mean_cov)[(48+length(significant_variable_lst)):(length(significant_variable_lst)+86)]
  
  colnames(low_CI_data) = colnames(MD_pheno_sig_volume_DKT_mean_cov)[48:(47+length(significant_variable_lst))]
  rownames(low_CI_data) = colnames(MD_pheno_sig_volume_DKT_mean_cov)[(48+length(significant_variable_lst)):(length(significant_variable_lst)+86)]
  colnames(up_CI_data) = colnames(MD_pheno_sig_volume_DKT_mean_cov)[48:(47+length(significant_variable_lst))]
  rownames(up_CI_data) = colnames(MD_pheno_sig_volume_DKT_mean_cov)[(48+length(significant_variable_lst)):(length(significant_variable_lst)+86)]
  rownames(betas_data)[2] = "Thalamus"
  rownames(betas_data)[8] = "Accumbens"
  rownames(pvalues_data)[2] = "Thalamus"
  rownames(pvalues_data)[8] = "Accumbens"
  
  rownames(low_CI_data)[2] = "Thalamus"
  rownames(low_CI_data)[8] = "Accumbens"
  rownames(up_CI_data)[2] = "Thalamus"
  rownames(up_CI_data)[8] = "Accumbens"
  
  rownames(betas_data) = capitalize(rownames(betas_data))
  rownames(pvalues_data) = capitalize(rownames(pvalues_data))
  
  #对FieldID进行替换成description Serum NMR
  
  tmp = colnames(betas_data)
  for (i in c(1:length(tmp))) {
    if(tmp[i] %in% NMR_Field$FieldID){
      tmp[i] = NMR_Field$Abbreviation[NMR_Field$FieldID == tmp[i]]
    }
    if(tmp[i] %in% Serum_Field$FieldID){
      tmp[i] = Serum_Field$Abbreviation[Serum_Field$FieldID == tmp[i]]
    }
  }
  colnames(betas_data) = tmp
  colnames(pvalues_data) = tmp
  colnames(low_CI_data) = tmp
  colnames(up_CI_data) = tmp
  
  pvalues_data_bfi = pvalues_data
  for (i in c(1:39)) {
    pvalues_data_bfi[i,] = p.adjust(pvalues_data[i,],
                                    method = "bonferroni")
  }
  
  #把全都不显著的列去掉
  b = c()
  for (i in (1:length(colnames(pvalues_data_bfi)))) {
    a = 0
    for (j in (1:length(rownames(pvalues_data_bfi)))) {
      if(pvalues_data_bfi[j,i] < 0.05){
        a = a+1
      }
    }
    if(a!=0){b=c(b,i)}
  }
  pvalues_data_bfi = as.data.frame(pvalues_data_bfi[,b])
  betas_data = as.data.frame(betas_data[,b])
  pvalues_all = cbind(pvalues_all, pvalues_data_bfi)
  betas_all = cbind(betas_all,betas_data)
  
  low_CI_data = as.data.frame(low_CI_data[,b])
  up_CI_data = as.data.frame(up_CI_data[,b])
  low_CI_all = cbind(low_CI_all, low_CI_data)
  up_CI_all = cbind(up_CI_all, up_CI_data)
  
}
betas_all = betas_all[-1]
pvalues_all = pvalues_all[-1]
low_CI_all = low_CI_all[-1]
up_CI_all =  up_CI_all[-1]

colnames(betas_all)[9] = "LYMP-P"
colnames(pvalues_all)[9] = "LYMP-P"
colnames(low_CI_all)[9] = "LYMP-P"
colnames(up_CI_all)[9] = "LYMP-P"
pvalues_matrix = as.matrix(pvalues_all)
betas_matrix = as.matrix(betas_all)

p <- melt(pvalues_matrix)
d <- melt(betas_matrix)

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
                    "*** p_bfi < 0.0005","\n\n","beta"))    

setwd("/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/Figures/brain_blood_corr/")
ggsave(paste("DKT_all_sig.pdf",sep = ""),width = 10,height = 10)
