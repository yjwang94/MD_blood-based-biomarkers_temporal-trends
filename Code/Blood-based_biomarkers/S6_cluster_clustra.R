#install.packages("clustra")
library(clustra)
library(reshape2)
loess_predict_mean_all = read.csv("/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/CaseControl/S5_loess_predict_mean_all.csv")
loess_predict_mean = loess_predict_mean_all[,32:80]
#AX:1:31. DEP:32:80. BIP:81:90. SCH:91:103
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
colnames(loess_predict_mean) = c("time","id","response")


#Fuzzy C-Means Clustering
#########AX:
#0.8-10
res.kmeans <- clustra(loess_predict_mean, k=25 ,maxdf = 10, conv = c(1000, 0))
res.kmeans$group
tabulate(res.kmeans$group)
sil = clustra_sil(res.kmeans)

lapply(sil, plot_silhouette)

set.seed(331081199809040020)
res.kmeans <- clustra(loess_predict_mean, k=18 ,maxdf = 5, conv = c(1000, 0))
res.kmeans$group
tabulate(res.kmeans$group)
sil = clustra_sil(res.kmeans)

lapply(sil, plot_silhouette)

#########DEP:
#0.84-14
res.kmeans <- clustra(loess_predict_mean, k=26 ,maxdf = 15, conv = c(1000, 0))
res.kmeans$group
tabulate(res.kmeans$group)
sil = clustra_sil(res.kmeans)

lapply(sil, plot_silhouette)

#0.83-13
set.seed(123)
res.kmeans <- clustra(loess_predict_mean, k=19 ,maxdf = 6, conv = c(1000, 0))
res.kmeans$group
tabulate(res.kmeans$group)
sil = clustra_sil(res.kmeans)

lapply(sil, plot_silhouette)


#########BIP:
#0.95-8
set.seed(1)
res.kmeans <- clustra(loess_predict_mean, k=19 ,maxdf = 6, conv = c(10000, 0))
res.kmeans$group
tabulate(res.kmeans$group)
sil = clustra_sil(res.kmeans)

lapply(sil, plot_silhouette)

#0.92-5
set.seed(12)
res.kmeans <- clustra(loess_predict_mean, k=10 ,maxdf = 6, conv = c(10000, 0))
res.kmeans$group
tabulate(res.kmeans$group)
sil = clustra_sil(res.kmeans)

lapply(sil, plot_silhouette)

#########SCH:
#0.92-10
set.seed(12345678)
res.kmeans <- clustra(loess_predict_mean, k=10 ,maxdf = 5, conv = c(1000, 0))
res.kmeans$group
tabulate(res.kmeans$group)
sil = clustra_sil(res.kmeans)

lapply(sil, plot_silhouette)
#0.92-10
set.seed(12345678)
res.kmeans <- clustra(loess_predict_mean, k=20 ,maxdf = 5, conv = c(1000, 0))
res.kmeans$group
tabulate(res.kmeans$group)
sil = clustra_sil(res.kmeans)

lapply(sil, plot_silhouette)
#0.92-10
set.seed(1)
res.kmeans <- clustra(loess_predict_mean, k=12 ,maxdf = 10, conv = c(1000, 0))
res.kmeans$group
tabulate(res.kmeans$group)
sil = clustra_sil(res.kmeans)

lapply(sil, plot_silhouette)

#########################
cluster_data = as.data.frame(res.kmeans$group)
rownames(cluster_data) = cmp
colnames(cluster_data) = "cluster"
tmp = rownames(cluster_data)
tmp <- substr(tmp, 1, 9)
rownames(cluster_data) = tmp
for (i in (1:length(rownames(cluster_data)))) {
  tmp[i] = data_FieldID$Abbreviation[data_FieldID$FieldID == tmp[i]]
}
rownames(cluster_data) = tmp
write.csv(cluster_data,"/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/CaseControl/S6_cluster_clustra/S62_cluster_data_SCH1.csv")

