#提取脑区数据
library(stringr)
library(reshape2)
library(dplyr)
library(ggplot2)
#install.packages("htmltools")
library(Hmisc)

MD = "DEP"

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

#############
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

#MD血液指标和脑区体积进行配对， with 39个脑区体积 +协变量
#并且剔除其他疾病的人群53903
MD_zscore = read.csv(paste("/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/CaseControl/",MD,"/S4_zdata_MD.csv",sep = ""))
disease_exclude_eid = read.csv("/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/Neuroimaging/disease exclude eid.csv")
disease_exclude_eid_lst = disease_exclude_eid$eid
MD_zscore = MD_zscore[!(MD_zscore$eid %in% disease_exclude_eid_lst),]
significant_variable = read.csv(paste("/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/CaseControl/S6_cluster_clustra/S62_cluster_data_",MD,".csv",sep = ""))
significant_variable_sorted <- arrange(significant_variable, cluster)
significant_variable_lst = significant_variable$X
significant_variable_sorted_lst = significant_variable_sorted$X
# rename pheno of UKB_pheno
colnames(MD_zscore)[48:(47+length(significant_variable_lst))] = significant_variable_lst
#对MD_zscore中的blood变量按照cluster排序
MD_zscore_arrange = cbind(MD_zscore[,1:47],MD_zscore[,significant_variable_sorted_lst])
#把脑区体积和变量数据进行匹配按照eid
MD_pheno_sig_volume_DKT_mean_cov = merge(MD_zscore_arrange, DTI_mean, on = "eid")


#拟合线性模型
betas_data = data.frame(matrix(nrow = 39))
pvalues_data = data.frame(matrix(nrow = 39))
for (i in (48:(47+length(significant_variable_lst)))){#blood biomarkers
  betas = data.frame()
  pvalues = data.frame()
  for (j in ((48+length(significant_variable_lst)):(length(significant_variable_lst)+86))) {#脑区
    complete = MD_pheno_sig_volume_DKT_mean_cov[complete.cases(MD_pheno_sig_volume_DKT_mean_cov[,i]),]
    lm1 = lm(complete[,i] ~ complete[,'group']*complete[,j]+ Gender + TD_Index + Age + SMK_Status1 + SMK_Status2  +
               ALC_Status1 + ALC_Status2 + Education3 + Education4 + 
               Education5 + Education7 + Education8 + 
               Ethnicity1 + Ethnicity2 + Ethnicity3 + Sites11001 + Sites11002 + Sites11003 +
               Sites11004 + Sites11005 + Sites11006 + Sites11007 + Sites11008 + Sites11009 +
               Sites11010 + Sites11011 + Sites11012 + Sites11013 + Sites11014 + Sites11016 + 
               Sites11017 + Sites11018 + Sites11020 + Sites11021 + Sites11022 + Sites11023, data = complete)
    beta = lm1$coefficients[2]
    sumres <- summary(lm1)
    pvalue = sumres$coefficients[, "Pr(>|t|)"][2]
    betas = rbind(betas,beta)
    pvalues = rbind(pvalues,pvalue)
  }
  betas_data = cbind(betas_data, betas)
  pvalues_data = cbind(pvalues_data, pvalues)
}

betas_data = betas_data[,-1]
pvalues_data = pvalues_data[,-1]
colnames(betas_data) = colnames(MD_pheno_sig_volume_DKT_mean_cov)[48:(47+length(significant_variable_lst))]
rownames(betas_data) = colnames(MD_pheno_sig_volume_DKT_mean_cov)[(48+length(significant_variable_lst)):(length(significant_variable_lst)+86)]
colnames(pvalues_data) = colnames(MD_pheno_sig_volume_DKT_mean_cov)[48:(47+length(significant_variable_lst))]
rownames(pvalues_data) = colnames(MD_pheno_sig_volume_DKT_mean_cov)[(48+length(significant_variable_lst)):(length(significant_variable_lst)+86)]
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

pvalues_matrix = as.matrix(pvalues_data_bfi)
betas_matrix = as.matrix(betas_data)

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
ggsave(paste("DTI_",MD,"_sig.pdf",sep = ""),width = 6,height = 8)
