library(ggplot2)
library(reshape2)
library(dplyr)
library(cluster)
library(NbClust)

MD = "DEP"
fpath = paste('/data/Mental_assessment/', MD, '/S1_case_control_match_yrs_online_date.csv', sep = '')
data = read.csv(fpath)
colnames(data)[2] = "duration_FL"
covariance = read.csv("/data/raw_data/Covariates_Data_Imputed.csv")
MD_score = read.csv("/data/raw_data/UKB_mental.csv")

#zdata是下面的代码准备的，为了不再重新跑一遍直接进行存储
data_cov = merge(data, covariance, by = 'eid')


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

data_cov_pheno = merge(data_cov, MD_score, by = 'eid')


#eliminate the influence of covariance
data_cov_pheno_elimate = data_cov_pheno
for (i in (48:length(colnames(data_cov_pheno)))){
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
for (i in (48:length(colnames(data_cov_pheno)))) {
  zdata[complete.cases(zdata[,i]),][,i] = (zdata[complete.cases(zdata[,i]),][,i] - mean(zdata[complete.cases(zdata[,i]),][,i]))/sd(zdata[complete.cases(zdata[,i]),][,i])
}
#colnames(zdata)[2] = 'duration'

write.csv(zdata, paste("/data/Mental_assessment/", MD, '/S2_zdata_mental.csv',sep = ''), row.names = FALSE)


