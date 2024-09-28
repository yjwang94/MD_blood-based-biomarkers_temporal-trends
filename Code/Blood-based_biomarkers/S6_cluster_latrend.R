#install.packages("latrend")
#install.packages("qqplotr")
library(latrend)
library(reshape2)
loess_predict_mean_all = read.csv("/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/CaseControl/S5_loess_predict_mean_all.csv")
loess_predict_mean = loess_predict_mean_all[,42:90]
# rename pheno of UKB_pheno
tmp = colnames(loess_predict_mean)[1:length(colnames(loess_predict_mean))]
for (i in (1:length(colnames(loess_predict_mean)))) {
  tmp[i] = sub('X', '', tmp[i])
  tmp[i] = sub('\\.', '\\-', tmp[i])
}
colnames(loess_predict_mean)[1:length(colnames(loess_predict_mean))] = tmp

data_FieldID = read.csv("/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/PhenoData/data_FieldID_abbre.csv")
tmp = colnames(loess_predict_mean)[1:length(colnames(loess_predict_mean))]
dmp = tmp
cmp = tmp
hmp = tmp
for (i in (1:length(colnames(loess_predict_mean)))) {
  tmp[i] = substr(tmp[i], 1, 9)
  dmp[i] = substring(cmp[i], 10)
  tmp[i] = data_FieldID$Abbreviation[data_FieldID$FieldID == tmp[i]]
  hmp[i] = paste(tmp[i],dmp[i], sep = "_")
}
colnames(loess_predict_mean)[1:length(colnames(loess_predict_mean))] = hmp
loess_predict_mean$time = c(-10:9)
#变成长矩阵
loess_predict_mean <- melt(loess_predict_mean,id.vars = c("time"))
colnames(loess_predict_mean) = c("Time","Id","Y")

method <- lcMethodLMKM(Y ~ Time, id = "Id", time = "Time")
model <- latrend(method, data = loess_predict_mean, nClusters = 12, seed = 1)
plot(model)
cluster_data = as.data.frame(model@model$cluster)
rownames(cluster_data) = hmp
colnames(cluster_data) = "cluster"
tmp = rownames(cluster_data)
tmp <- sub("_[A-Z]+$", "", tmp)
rownames(cluster_data) = tmp
for (i in (1:length(rownames(cluster_data)))) {
  tmp[i] = data_FieldID$FieldID[data_FieldID$Abbreviation == tmp[i]]
}
rownames(cluster_data) = tmp
write.csv(cluster_data,"/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/CaseControl/S6_cluster_latrend/S62_cluster_data_DEP12.csv")
