#DTI
library(stringr)
fpath = "/data/Neuroimaging/"
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

#MD血液指标和脑区体积进行配对，所有MD（no control）显著指标的并集 with 48个变量
#所有MD并且剔除其他疾病的人群53903
MD_lst = c("AX", "BIP", "DEP", "SCH")
MD_case = data.frame()
for (MD in MD_lst) {
  case_path = paste("/data/CaseControl/",MD,"/S1_",MD,"_case_control_eid_df.csv",sep = "")
  MD_data = read.csv(case_path)
  MD_case = rbind(MD_case, MD_data)
}
MD_case_lst = unique(MD_case$case_id)
disease_exclude_eid = read.csv("/data/Neuroimaging/disease exclude eid.csv")
disease_exclude_eid_lst = disease_exclude_eid$eid
MD_case_lst = setdiff(MD_case_lst, disease_exclude_eid_lst)

#匹配Serum & NMR
UKB_pheno = read.csv("/data/PhenoData2Analysis/UKB_pheno.csv")
UKB_blood_dict = read.csv("/data/PhenoData2Analysis/UKB_Blood_dict.csv")
NMR_dict = read.csv("/data/PhenoData2Analysis/abbreviation_NMR.csv")
Serum_dict = read.csv("/data/PhenoData2Analysis/Serum_dict1.csv")
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
#所有MD个体Serum和NMR的值
MD_pheno = UKB_pheno227[UKB_pheno227$eid %in% MD_case_lst,]

#MD显著变量的并集
MD_lst = c("AX", "BIP", "DEP", "SCH")
significant_variable = data.frame()
for (MD in MD_lst) {
  variable_path = paste("/data/CaseControl/S6_cluster_data_",MD,".csv",sep = "")
  variable_data = read.csv(variable_path)
  significant_variable = rbind(significant_variable, variable_data)
}
significant_variable_lst = significant_variable$X
for (i in c(1:95)) {
  significant_variable_lst[i] = substr(significant_variable_lst[i],1,9)
}
significant_variable_lst = unique(significant_variable_lst)

MD_pheno_sig = cbind(MD_pheno[1],MD_pheno[colnames(MD_pheno) %in% significant_variable_lst])

#把脑区和变量数据进行匹配按照eid
MD_pheno_sig_volume_DKT_mean = merge(MD_pheno_sig, DTI_mean, on = "eid")

#对脑区体积和变量数据剔除协变量处理 对连续变量进行规范化不太影响就不做了
cov = read.csv("/data/CaseControl/Covariates_Data_Imputed.csv")
MD_pheno_sig_volume_DKT_mean_cov = merge(MD_pheno_sig_volume_DKT_mean, cov, by = 'eid')
#哑变量处理
MD_pheno_sig_volume_DKT_mean_cov$SMK_Status = as.character(MD_pheno_sig_volume_DKT_mean_cov$SMK_Status)
MD_pheno_sig_volume_DKT_mean_cov$ALC_Status = as.character(MD_pheno_sig_volume_DKT_mean_cov$ALC_Status)
MD_pheno_sig_volume_DKT_mean_cov$Gender = as.character(MD_pheno_sig_volume_DKT_mean_cov$Gender)
MD_pheno_sig_volume_DKT_mean_cov$Education = as.character(MD_pheno_sig_volume_DKT_mean_cov$Education)
MD_pheno_sig_volume_DKT_mean_cov$Ethnicity = as.character(MD_pheno_sig_volume_DKT_mean_cov$Ethnicity)
MD_pheno_sig_volume_DKT_mean_cov$Sites = as.character(MD_pheno_sig_volume_DKT_mean_cov$Sites)

dumy <- model.matrix( ~ SMK_Status - 1, data = MD_pheno_sig_volume_DKT_mean_cov)
MD_pheno_sig_volume_DKT_mean_cov <- cbind(MD_pheno_sig_volume_DKT_mean_cov,dumy[,-1])
dumy <- model.matrix( ~ ALC_Status - 1, data = MD_pheno_sig_volume_DKT_mean_cov)
MD_pheno_sig_volume_DKT_mean_cov <- cbind(MD_pheno_sig_volume_DKT_mean_cov,dumy[,-1])
dumy <- model.matrix( ~ Education - 1, data = MD_pheno_sig_volume_DKT_mean_cov)
MD_pheno_sig_volume_DKT_mean_cov <- cbind(MD_pheno_sig_volume_DKT_mean_cov,dumy[,-1])
dumy <- model.matrix( ~ Ethnicity - 1, data = MD_pheno_sig_volume_DKT_mean_cov)
MD_pheno_sig_volume_DKT_mean_cov <- cbind(MD_pheno_sig_volume_DKT_mean_cov,dumy[,-1])
dumy <- model.matrix( ~ Sites - 1, data = MD_pheno_sig_volume_DKT_mean_cov)
MD_pheno_sig_volume_DKT_mean_cov <- cbind(MD_pheno_sig_volume_DKT_mean_cov,dumy[,-1])

MD_pheno_sig_volume_DKT_mean_cov_elimate = MD_pheno_sig_volume_DKT_mean_cov
for (i in (2:95)){
  complete = MD_pheno_sig_volume_DKT_mean_cov_elimate[complete.cases(MD_pheno_sig_volume_DKT_mean_cov_elimate[,i]),]
  glm1 = glm(complete[,i] ~ Gender + TD_Index + Age + SMK_Status1 + SMK_Status2 +
               ALC_Status1 + ALC_Status2 + Education3 + Education4 + 
               Education5 + Education7 + Education8 + 
               Ethnicity1 + Ethnicity2 + Ethnicity3 + Sites11001 + Sites11002 + Sites11003 +
               Sites11004 + Sites11005 + Sites11006 + Sites11007 + Sites11008 + Sites11009 + 
               Sites11010 + Sites11011 + Sites11012 + Sites11013 + Sites11014 + Sites11016 + 
               Sites11017 + Sites11018 + Sites11020 + Sites11021 + Sites11022 + Sites11023, data = complete)
  
  res = residuals(glm1)
  
  new_data = res+glm1$coefficients[1]
  
  
  MD_pheno_sig_volume_DKT_mean_cov_elimate[complete.cases(MD_pheno_sig_volume_DKT_mean_cov_elimate[,i]),][,i] = new_data
} 

#画热图
library(reshape2)
library(dplyr)
library(ggplot2)

corrs_data = data.frame(matrix(nrow = 42))
pvalues_data = data.frame(matrix(nrow = 42))
for (i in c(2:53)) {
  corrs = data.frame()
  pvalues = data.frame()
  for (j in c(54:95)) {
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
colnames(corrs_data) = colnames(MD_pheno_sig_volume_DKT_mean_cov_elimate)[2:53]
rownames(corrs_data) = colnames(MD_pheno_sig_volume_DKT_mean_cov_elimate)[54:95]
colnames(pvalues_data) = colnames(MD_pheno_sig_volume_DKT_mean_cov_elimate)[2:53]
rownames(pvalues_data) = colnames(MD_pheno_sig_volume_DKT_mean_cov_elimate)[54:95]
#rownames(corrs_data)[2] = "Thalamus"
#rownames(corrs_data)[8] = "Accumbens"
library(Hmisc)
#rownames(corrs_data) = capitalize(rownames(corrs_data))
#rownames(pvalues_data) = capitalize(rownames(pvalues_data))

#对FieldID进行替换成description Serum NMR
tmp = colnames(corrs_data)
for (i in c(1:52)) {
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
    0.05 < value1 & value1 < 0.5 ~ paste("*"), 
    0.005 < value1 & value1 <= 0.05 ~ paste("**"),
    value1 <= 0.005 ~ paste("***")))
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

