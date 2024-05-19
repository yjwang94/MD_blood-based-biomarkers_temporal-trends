library(ggplot2)
library(reshape2)
library(dplyr)
library(cluster)
library(NbClust)
library(RColorBrewer)


MD = c("Anxiety", "Bipolar", "Depression", "Schizophrenia")
i = "Schizophrenia"
zdata_all = data.frame()
for (i in MD) {
  zdata = read.csv(paste("/data/Mental_assessment/", i, '/S2_zdata_mental.csv',sep = ''))
  zdata$MD = i
  zdata1 = zdata[c("duration_FL","anxiety",	"mania",	"wellbeing",	"psychotic_experience",	"selfharm",	
                   "mental_distress",	"trauma",	"depressnew","MD")]
  zdata_all = rbind(zdata_all, zdata1)
}
zdata_all$duration_FL = zdata_all$duration_FL*(-1)

setwd("/data/Mental_assessment/")
color = c('darkred', 'chocolate3', 'chartreuse3', 
          'deepskyblue3', 'purple3','cyan4',
          'navajowhite2', 'navajowhite4',
          'seagreen3','burlywood3','coral3','violetred2',
          'cyan3','magenta3','maroon','mediumorchid4','cadetblue2')
MD_score_lst = c("anxiety",	"mania",	"wellbeing",	"psychotic_experience",	"selfharm",	
                 "mental_distress",	"trauma",	"depressnew")

for (i in MD_score_lst) {
  zdata_plot_cluster = zdata[,c("duration_FL",i,"MD")]
  zdata_plot_melt = zdata_all[,c("duration_FL",i,"MD")]
  zdata_plot_melt = zdata_plot_melt[complete.cases(zdata_plot_melt[,2]),]
  colnames(zdata_plot_melt)[2] = "value"
  p <- ggplot(data=zdata_plot_melt, aes(x=duration_FL, y=value, group = MD)) +
    #geom_line(stat="smooth", data=zdata_plot_melt, aes(x=duration, y=i, group=MD),method = 'loess', span = 2,se = FALSE,size = 0.4, color = color[i],alpha = 0.4) +
    labs(x="Years to Diagnosis", y=i) +
    scale_x_continuous(limits=c(-10,10),breaks = seq(-10, 10, by = 2)) +
    # xy轴标签
    #geom_point() +
    # 添加散点
    geom_smooth(aes(color = MD),method = 'loess', span = 2, se = FALSE,linewidth = 6, show.legend = TRUE)+
    # 添加线性回归直线
    #theme(panel.grid=element_blank(), panel.background=element_rect(fill='transparent', color='black')) +
    # 去点网格、去掉背景、添加边框
    guides(color=guide_legend(title='Field'))+
    # 修改legend标题
    #scale_color_manual(values = c('darkred'))+
    #geom_hline(yintercept= 1, color = "#A9A9A9", lwd = 1.6)+
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
          axis.line.y=element_line(size=2)) + 
    scale_color_brewer(palette = "Set1")
  #legend.position="none")
  #修改配色
  #p_cluster
  ggsave(paste(i , ".pdf",sep = ""))
}

