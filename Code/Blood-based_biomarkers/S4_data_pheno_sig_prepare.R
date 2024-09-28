covariance = read.csv("/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/raw_data/Covariates_Data_Imputed.csv")
UKB_pheno = read.csv("/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/PhenoData2Analysis/UKB_pheno.csv")

MD = "DEP"
fpath = paste('/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/CaseControl/', MD, '/S1_case_control_match_yrs.csv', sep = '')
data = read.csv(fpath)
colnames(data)[2] = "duration"
pvalue = read.csv(paste("/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/CaseControl/",MD,"/S3_Pvalue_beta.csv",sep = ""))
pvalue_FiledID = pvalue$FieldID_full[pvalue$p_group_bfi<0.01 | pvalue$p_group_time_bfi< 0.01 | pvalue$p_group_time2_bfi< 0.01]
UKB_pheno_sig = UKB_pheno[, c('eid', pvalue_FiledID)]
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
continues = c("TD_Index", "Age", "BMI")

#data_cov1 = data_cov
for (i in continues) {
  data_cov[,i] = (data_cov[,i]-mean(data_cov[,i])) / sd(data_cov[,i])
}

UKB_pheno_sig = UKB_pheno[, c('eid', pvalue_FiledID)]

data_pheno_sig = merge(data_cov, UKB_pheno_sig, by = "eid")
#标准化
data_pheno_sig1 = as.data.frame(data_pheno_sig[48:length(colnames(data_pheno_sig))])
data_pheno_sig2 = cbind(data_pheno_sig[1:47],data_pheno_sig1)
write.csv(data_pheno_sig2,paste("/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/CaseControl/",MD,"/S4_data_pheno_sig.csv",sep = ""),row.names = FALSE)
