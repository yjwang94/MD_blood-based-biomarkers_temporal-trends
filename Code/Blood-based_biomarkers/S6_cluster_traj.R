library(traj)
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
loess_predict_mean = t(loess_predict_mean)
loess_predict_mean = as.data.frame(loess_predict_mean)
loess_predict_mean$ID = rownames(loess_predict_mean)
loess_predict_mean <- cbind(loess_predict_mean$ID, loess_predict_mean[, 1:20])
colnames(loess_predict_mean)[1] = "ID"

m = Step1Measures(loess_predict_mean, ID = TRUE, measures = 1:18, midpoint = 3)
s = Step2Selection(m)
s2 = Step2Selection(m)

c.part <- Step3Clusters(s2, nclusters = 10)$partition
cluster_data = c.part
colnames(cluster_data) = c("X","cluster")
tmp = cluster_data$X
tmp <- sub("_[A-Z]+$", "", tmp)
rownames(cluster_data) = tmp
for (i in (1:length(rownames(cluster_data)))) {
  tmp[i] = data_FieldID$FieldID[data_FieldID$Abbreviation == tmp[i]]
}
cluster_data$X = tmp
rownames(cluster_data) = tmp
write.csv(cluster_data,"/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/CaseControl/S6_cluster_traj/S62_cluster_data_DEP10.csv",row.names=FALSE)
