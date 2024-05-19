library(ggplot2)
library(reshape2)
library(dplyr)
library(cluster)
library(NbClust)

MD = "DEP"
UKB_pheno = read.csv("/data/PhenoData2Analysis/UKB_pheno.csv")
UKB_blood_dict = read.csv("/data/PhenoData2Analysis/UKB_Blood_dict.csv")
NMR_dict = read.csv("/data/PhenoData2Analysis/abbreviation_NMR.csv")
Serum_dict = read.csv("/data/PhenoData2Analysis/Serum_dict1.csv")

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


#对每个变量的loess的predict的数字进行聚类
zdata = read.csv(paste("/data/CaseControl/", MD, '/S4_zdata_MD.csv',sep = ''))
tmp = colnames(zdata)[48:length(colnames(zdata))]
for (i in (1:length(colnames(zdata)))) {
  tmp[i] = sub('X', '', tmp[i])
  tmp[i] = sub('\\.', '\\-', tmp[i])
}
colnames(zdata)[48:length(colnames(zdata))] = tmp
tmp = colnames(zdata)[48:length(colnames(zdata))]
zdata$duration =zdata$duration*(-1)

clustered_data = read.csv(paste("/data/CaseControl/S6_cluster_data_",MD,"_33.csv",sep = ""))

tmp = clustered_data$X
for (i in (1:length(colnames(zdata)))) {
  tmp[i] = sub('X', '', tmp[i])
  tmp[i] = sub('\\.', '\\-', tmp[i])
}
tmp = substring(tmp,1,9)
tmp = tmp[1:length(clustered_data$X)]
clustered_data$X = tmp


setwd(paste("/data/CaseControl/",MD,"/trend_cluster/",sep = ''))
color = c('darkred', 'chocolate3', 'chartreuse3', 
          'deepskyblue3', 'purple3','cyan4',
          'navajowhite2', 'navajowhite4',
          'seagreen3','burlywood3','coral3','violetred2',
          'cyan3','magenta3','maroon','mediumorchid4','cadetblue2')
i=2
for (i in c(1,2)) {
  zdata_plot_cluster = zdata[colnames(zdata) %in% 
                               c("duration",clustered_data[clustered_data$Cluster == i,]$X)]
  #zdata_plot_cluster2 = zdata[colnames(zdata) %in% 
                               #c(clustered_data[clustered_data$Cluster == 5,]$X)]
  #zdata_plot_cluster = cbind(zdata_plot_cluster1,zdata_plot_cluster2)
  tmp = colnames(zdata_plot_cluster)[2:length(colnames(zdata_plot_cluster))]
  for (m in (1:length(colnames(zdata_plot_cluster))-1)) {
    j = data_FieldID$Abbreviation[data_FieldID$FieldID == tmp[m]]
    tmp[m] = j
  }
  colnames(zdata_plot_cluster)[2:length(colnames(zdata_plot_cluster))] = tmp
  zdata_plot_melt = melt(zdata_plot_cluster, id = 'duration')
  zdata_plot_melt = zdata_plot_melt[complete.cases(zdata_plot_melt[,3]),]
  
  p_cluster <- ggplot() +
    geom_line(stat="smooth", data=zdata_plot_melt, aes(x=duration, y=value, group=variable),method = 'loess', span = 2,se = FALSE,linewidth = 1, color = color[i+4],alpha = 0.4) +
    labs(x="Years to Diagnosis", y="Z-Score") +
    scale_x_continuous(limits=c(-10,10),breaks = seq(-10, 10, by = 2)) +
    geom_smooth(data=zdata_plot_melt, aes(x=duration, y=value),method = 'loess', span = 2, se = FALSE,linewidth = 6, show.legend = TRUE,color = color[i+4])+
    #theme(panel.grid=element_blank(), panel.background=element_rect(fill='transparent', color='black')) +
    # 去点网格、去掉背景、添加边框
    guides(color=guide_legend(title='Field'))+
    scale_color_manual(values = c('darkred'))+
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
  ggsave(paste('cluster',i , ".pdf",sep = ""))
}


