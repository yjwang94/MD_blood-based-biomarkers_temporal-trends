library(ggplot2)
library(reshape2)
library(dplyr)
library(cluster)
library(NbClust)

MD = "AX"
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


raw_data = read.csv(paste("/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/CaseControl/",MD,"/S2_data_cov_pheno_forPvalue.csv",sep = ''))
tmp = colnames(raw_data)[15:length(colnames(raw_data))]
for (i in (1:length(colnames(raw_data)))) {
  tmp[i] = sub('X', '', tmp[i])
  tmp[i] = sub('\\.', '\\-', tmp[i])
  tmp[i] = substr(tmp[i],1,9)
}
colnames(raw_data)[15:length(colnames(raw_data))] = tmp
tmp = colnames(raw_data)[15:length(colnames(raw_data))]
for (m in (1:length(colnames(raw_data)))) {
  j = data_FieldID$Abbreviation[data_FieldID$FieldID == tmp[m]]
  tmp[m] = j
}
colnames(raw_data)[15:length(colnames(raw_data))] = tmp
raw_data$duration = raw_data$duration*(-1)
setwd(paste("//Users/yujia/Desktop/a/project/Evolution/MD/MD_review/Figures/raw_data_distribution/",MD,"/",sep = ''))

i=37

for (i in c(15:length(colnames(raw_data)))) {
  raw_data1 = raw_data[complete.cases(raw_data[,i]),]
  hist(raw_data1[,i]) 
    #geom_histogram()+
    #coord_cartesian(ylim=c(0,10000))
    #labs(x="Years to Diagnosis", y=colnames(raw_data1)[i])
    #scale_x_continuous(limits=c(-10,10),breaks = seq(-10, 10, by = 2)) +
    #geom_smooth(data=zdata_plot_melt, aes(x=duration, y=Field),method = 'loess', span = 2, se = FALSE,linewidth = 6, show.legend = TRUE,color = 'darkred')+
    #theme(panel.grid=element_blank(), panel.background=element_rect(fill='transparent', color='black')) +
    # 去点网格、去掉背景、添加边框
    #guides(color=guide_legend(title='Field'))+
    #scale_color_manual(values = c('darkred'))+
  #legend.position="none")
  #修改配色
  #p_cluster
  png(paste(colnames(raw_data)[i],".png",sep = ""),height=8,width=8)
  dev.off()
}

#把病人按病程分几组，画boxplot
#分成五组，每组四年
i=37
raw_data1 = raw_data[complete.cases(raw_data[,i]),]
raw_data1$class[raw_data1$duration >= (-10) & raw_data1$duration<(-6)] = "-10~-6"
raw_data1$class[raw_data1$duration >= (-6) & raw_data1$duration<(-2)] = "-6~-2"
raw_data1$class[raw_data1$duration >= (-2) & raw_data1$duration<(2)] = "-2~2"
raw_data1$class[raw_data1$duration >= (2) & raw_data1$duration<(6)] = "2~6"
raw_data1$class[raw_data1$duration >= (6) & raw_data1$duration<=(10)] = "6~10"

raw_data1$class = as.factor(raw_data1$class)
data_WBC = raw_data1[raw_data1$WBC<15,]
boxplot(data_WBC[,183] ~ class, #代表以count为x轴，以spray为y轴
        data = data_WBC,  #指定数据集为data
        col = "pink")
boxplot(raw_data1[,183] ~ class, #代表以count为x轴，以spray为y轴
        data = raw_data1,  #指定数据集为data
        col = "pink")
boxplot(raw_data1[,19] ~ class, #代表以count为x轴，以spray为y轴
        data = raw_data1,  #指定数据集为data
        col = "pink")
