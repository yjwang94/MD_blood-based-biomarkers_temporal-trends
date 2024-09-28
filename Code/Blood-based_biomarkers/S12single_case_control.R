library(ggplot2)
library(reshape2)
library(dplyr)
library(cluster)
library(NbClust)

MD = "BIP"
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


zdata = read.csv(paste("/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/CaseControl/", MD, '/S4_zdata_MD.csv',sep = ''))
tmp = colnames(zdata)[48:length(colnames(zdata))]
for (i in (1:length(colnames(zdata)))) {
  tmp[i] = sub('X', '', tmp[i])
  tmp[i] = sub('\\.', '\\-', tmp[i])
  tmp[i] = substr(tmp[i],1,9)
}
colnames(zdata)[48:length(colnames(zdata))] = tmp
tmp = colnames(zdata)[48:length(colnames(zdata))]
zdata$duration =zdata$duration*(-1)


setwd(paste("/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/Figures/single_case_control/",MD,"/",sep = ''))
color = c('darkred', 'chocolate3', 'chartreuse3', 
          'deepskyblue3', 'purple3','cyan4',
          'navajowhite2', 'navajowhite4',
          'seagreen3','burlywood3','coral3','violetred2',
          'cyan3','magenta3','maroon','mediumorchid4','cadetblue2')

for (i in c(48:length(colnames(zdata)))) {
  zdata_plot_cluster = zdata[c("duration","group",colnames(zdata)[i])]
  zdata_plot_cluster$duration = zdata_plot_cluster$duration*(-1)
  zdata_plot_cluster$group = as.factor(zdata_plot_cluster$group)
  tmp = colnames(zdata_plot_cluster)[3:length(colnames(zdata_plot_cluster))]
  for (m in (1:length(colnames(zdata_plot_cluster))-2)) {
    j = data_FieldID$Abbreviation[data_FieldID$FieldID == tmp[m]]
    tmp[m] = j
  }
  colnames(zdata_plot_cluster)[3:length(colnames(zdata_plot_cluster))] = tmp
  zdata_plot_melt = zdata_plot_cluster[complete.cases(zdata_plot_cluster[,3]),]
  colnames(zdata_plot_melt)[3] = "Field"
  p_cluster <- ggplot() +
    geom_line(stat="smooth", data=zdata_plot_melt, aes(x=duration, y=Field, group=group,color = group),method = 'loess', span = 2,linewidth = 6,se = FALSE) +
    labs(x="Years to Diagnosis", y=colnames(zdata_plot_cluster)[3]) +
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
  ggsave(paste(colnames(zdata_plot_cluster)[3], ".pdf",sep = ""))
}


