library(ggplot2)
library(reshape2)
library(dplyr)
library(cluster)
library(NbClust)

MD = "DEP"
covariance = read.csv("/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/raw_data/Covariates_Data_Imputed.csv")
UKB_pheno = read.csv("/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/PhenoData2Analysis/UKB_pheno.csv")
UKB_blood_dict = read.csv("/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/PhenoData2Analysis/UKB_Blood_dict.csv")
NMR_dict = read.csv("/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/PhenoData2Analysis/abbreviation_NMR.csv")
Serum_dict = read.csv("/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/PhenoData2Analysis/Serum_dict1.csv")

# rename pheno of UKB_blood_dict
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

#读取数据并消除协变量的影响
fpath = paste('/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/CaseControl/', MD, '/S1_case_control_match_yrs.csv', sep = '')
data = read.csv(fpath)
colnames(data)[2] = "duration"
pvalue = read.csv(paste("/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/CaseControl/",MD,"/S3_Pvalue_beta.csv",sep = ""))
pvalue_FiledID = pvalue$FieldID_full[pvalue$p_group_bfi<0.01 | pvalue$p_group_time_bfi< 0.01 | pvalue$p_group_time2_bfi< 0.01]
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

data_cov_pheno = merge(data_cov, UKB_pheno_sig, by = 'eid')


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
raw_data = data_cov_pheno_elimate

tmp = colnames(raw_data)[48:length(colnames(raw_data))]
#for (i in (1:length(colnames(raw_data)))) {
#  tmp[i] = sub('X', '', tmp[i])
#  tmp[i] = sub('\\.', '\\-', tmp[i])
#  tmp[i] = substr(tmp[i],1,9)
#}
#colnames(raw_data)[15:length(colnames(raw_data))] = tmp
for (m in (1:length(colnames(raw_data)))) {
  j = data_FieldID$Abbreviation[data_FieldID$FieldID == tmp[m]]
  tmp[m] = j
}
colnames(raw_data)[48:length(colnames(raw_data))] = tmp
raw_data$duration = raw_data$duration*(-1)

setwd(paste("/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/Figures/single_case_control/",MD,"/",sep = ''))

for (i in c(48:length(colnames(raw_data)))) {
  raw_data1 = raw_data[c("duration","group",colnames(raw_data)[i])]
  raw_data1$group = as.factor(raw_data1$group)
  raw_data1 = raw_data1[complete.cases(raw_data1[,3]),]
  raw_data2 = raw_data1
  colnames(raw_data1)[3] = "Field"
  #raw_data1 = raw_data1[raw_data1$Field<150,]
  p_cluster <- ggplot(data=raw_data1,
                      aes(x=duration, y=Field, group=group,color = group)) + 
    geom_smooth(method = 'loess', span = 2, se = FALSE,size = 5,alpha=0.22,show.legend = TRUE)+
    labs(x="Years to Diagnosis", y=colnames(raw_data2)[3]) +
    scale_x_continuous(limits=c(-10,10),breaks = seq(-10, 10, by = 2)) +
    #geom_smooth(data=zdata_plot_melt, aes(x=duration, y=Field),method = 'loess', span = 2, se = FALSE,linewidth = 6, show.legend = TRUE,color = 'darkred')+
    #theme(panel.grid=element_blank(), panel.background=element_rect(fill='transparent', color='black')) +
    # 去点网格、去掉背景、添加边框
    #guides(color=guide_legend(title='Field'))+
    #scale_color_manual(values = c('darkred'))+
    scale_color_brewer(palette = "Set1")+
    theme_classic()+
    theme(plot.margin=unit(c(0.2,0.2,0.2,0.2), "cm"),
          panel.spacing =unit(c(0,0,0,0), "cm"),
          axis.text.x = element_text(color="black", size=40,face = "bold"),
          axis.text.y = element_text(color="black", size=40,face = "bold"),
          axis.title.x = element_text(color="black", size=40,face = "bold"),
          axis.title.y = element_text(color="black", size=40,face = "bold"),
          axis.ticks.x = element_line(colour = "black", size = 1.2),
          axis.ticks.y = element_line(colour = "black", size = 1.2),
          axis.ticks.length.x = unit(0.4,'cm'),
          axis.ticks.length.y = unit(0.4,'cm'),
          axis.line.x=element_line(size=2),
          axis.line.y=element_line(size=2))
  
  #
  #legend.position="none")
  #修改配色
  #p_cluster
  ggsave(paste(colnames(raw_data2)[3], ".pdf",sep = ""))
}



raw_data2$class[raw_data1$duration >= (-10) & raw_data1$duration<(-6)] = "-10~-6"
raw_data2$class[raw_data1$duration >= (-6) & raw_data1$duration<(-2)] = "-6~-2"
raw_data2$class[raw_data1$duration >= (-2) & raw_data1$duration<(2)] = "-2~2"
raw_data2$class[raw_data1$duration >= (2) & raw_data1$duration<(6)] = "2~6"
raw_data2$class[raw_data1$duration >= (6) & raw_data1$duration<=(10)] = "6~10"

boxplot(raw_data2[,3] ~ class, #代表以count为x轴，以spray为y轴
        data = raw_data2,  #指定数据集为data
        col = "pink")
library(parallel)
#Calculate the number of cores检查电脑当前可用核数
no_cores<-detectCores(logical=F)		#F-物理CPU核心数/T-逻辑CPU核心数

#Initiate cluster发起集群,同时创建数个R进行并行计算
#只是创建待用的核,而不是并行运算环境
cl<-makeCluster(no_cores)

#现只需要使用并行化版本的lapply,parLapply就可以
parLapply(cl, 1:10000,function(exponent) 2^exponent)
#当结束后要关闭集群,否则电脑内存会始终被R占用
stopCluster(cl)

