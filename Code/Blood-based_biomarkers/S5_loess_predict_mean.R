
MD_lst = c("AX","DEP","BIP","SCH")
loess_predict_mean_all = data.frame(matrix(nrow = 20))
MD = "BIP"
for (MD in MD_lst) {
  zdata1 = read.csv(paste("/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/CaseControl/", MD, '/S4_zdata_MD.csv',sep = ''))
  zdata1 = zdata1[zdata1$group == 1,]
  tmp = colnames(zdata1)[48:length(colnames(zdata1))]
  for (i in (1:length(colnames(zdata1)))) {
    tmp[i] = sub('X', '', tmp[i])
    tmp[i] = sub('\\.', '\\-', tmp[i])
  }
  colnames(zdata1)[48:length(colnames(zdata1))] = tmp
  tmp = colnames(zdata1)[48:length(colnames(zdata1))]
  
  loess_predict = list()
  zdata_forloess = zdata1[,c(2,48:length(colnames(zdata1)))]
  zdata_forloess$duration = zdata_forloess$duration*(-1)
  loess_predict_mean_MD = data.frame(matrix(nrow = 20))
  for(i in (2:length(zdata_forloess))){
    zdata_complete = zdata_forloess[complete.cases(zdata_forloess[,i]),]
    loess_variable <- loess(zdata_forloess[complete.cases(zdata_forloess[i]),][,i] ~ duration, data=zdata_forloess[complete.cases(zdata_forloess[i]),], span=2)
    smooth_variable <- predict(loess_variable)
    loess_predict = cbind(zdata_complete[1], as.data.frame(smooth_variable))
    #calculate mean loess predict value in 21 time frame
    loess_predict_mean = data.frame()
    for (j in (-10:9)) {
      m = mean(loess_predict[loess_predict$duration >= j & loess_predict$duration < j+1,][,2])
      loess_predict_mean = rbind(loess_predict_mean,m)
      rownames(loess_predict_mean)[j+11] = j
    }
    colnames(loess_predict_mean) = colnames(zdata_forloess)[i]
    loess_predict_mean_MD = cbind(loess_predict_mean_MD, loess_predict_mean)
  }
  loess_predict_mean_MD = loess_predict_mean_MD[-1]
  
  loess_predict_mean_all = cbind(loess_predict_mean_all, loess_predict_mean_MD)
}
loess_predict_mean_all = loess_predict_mean_all[-1]
#loess_predict_mean_AX = read.csv("/Users/yujia/Desktop/a/project/Evolution/MD/MD1205/S5_loess_predict_mean_AX.csv")
#loess_predict_mean_all = cbind(loess_predict_mean_AX,loess_predict_mean_all)
write.csv(loess_predict_mean_all, 
      "/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/CaseControl/S5_loess_predict_mean_all.csv",row.names = FALSE)

