loess_predict_mean_all = read.csv("/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/CaseControl/S5_loess_predict_mean_all.csv")
loess_predict_mean = loess_predict_mean_all[,32:80]
#AX:1:31. DEP:32:80. BIP:81:90. SCH:91:103
# rename pheno of UKB_pheno
tmp = colnames(loess_predict_mean)[1:length(colnames(loess_predict_mean))]
for (i in (1:length(colnames(loess_predict_mean)))) {
  tmp[i] = sub('X', '', tmp[i])
  tmp[i] = sub('\\.', '\\-', tmp[i])
  tmp[i] = substr(tmp[i], 1, 9)
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
  #hmp[i] = paste(tmp[i],dmp[i], sep = "_")
}
colnames(loess_predict_mean)[1:length(colnames(loess_predict_mean))] = tmp

#Hierarchical clustering dendrogram
loess_predict_mean = t(loess_predict_mean)
d <- dist(loess_predict_mean, method = "euclidean") 
fit2 <- hclust(d, method="ward.D") 
pdf("/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/CaseControl/S6_cluster_separate/S6_cluster_DEP2.pdf", width = 8, height = 6)
plot(fit2) 
groups <- cutree(fit2, k=10)
rect.hclust(fit2, k=10, border="red")
dev.off()

#The classification of the significant variables
clustered_data <- data.frame(DataPoint = 1:nrow(loess_predict_mean), Cluster = groups)
write.csv(clustered_data,"/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/CaseControl/S6_cluster_separate/S62_cluster_data_DEP2.csv")






cluster4 = read.csv("/data/CaseControl/S6_cluster_data_DEP.csv")

loess_predict_mean_DEP3 = loess_predict_mean_DEP[,colnames(loess_predict_mean_DEP) %in% cluster4$X[cluster4$Cluster == 3]]
loess_predict_mean = t(loess_predict_mean_DEP3)
d <- dist(loess_predict_mean, method = "euclidean") # 计算各个样本点之间的距离
fit2 <- hclust(d, method="ward.D") #进行Ward层次聚类
plot(fit2) # 绘制树状图展示聚类结果
groups <- cutree(fit2, k=3) # 设定聚类个数为3
# 给聚成的3个类别加上红色边框
rect.hclust(fit2, k=3, border="red")
clustered_data33 <- data.frame(DataPoint = 1:nrow(loess_predict_mean), Cluster = groups)
write.csv(clustered_data33,"/data/CaseControl/S6_cluster_data_DEP_33.csv")



cluster_data_AX = clustered_data[1:31,]
cluster_data_BIP = clustered_data[32:41,]
cluster_data_DEP = clustered_data[42:90,]
cluster_data_SCH = clustered_data[91:103,]
write.csv(cluster_data_AX, "/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/CaseControl/S63_cluster_data_AX.csv")
write.csv(cluster_data_BIP, "/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/CaseControl/S63_cluster_data_BIP.csv")
write.csv(cluster_data_DEP, "/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/CaseControl/S63_cluster_data_DEP.csv")
write.csv(cluster_data_SCH, "/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/CaseControl/S63_cluster_data_SCH.csv")
#write.csv(cluster_data_SLP, "/Users/yujia/Desktop/a/project/Evolution/MD/MD_prepost/zdata_trajectory/significant_variable/SLP/cluster_data_SLP.csv", row.names = FALSE)
#write.csv(cluster_data_SU, "/Users/yujia/Desktop/a/project/Evolution/MD/MD_prepost/zdata_trajectory/significant_variable/SU/cluster_data_SU.csv", row.names = FALSE)

