loess_predict_mean_all = read.csv("/data/CaseControl/S5_loess_predict_mean_all.csv")
loess_predict_mean_DEP = loess_predict_mean_all[,32:81]
loess_predict_mean = t(loess_predict_mean_DEP)
d <- dist(loess_predict_mean, method = "euclidean") # 计算各个样本点之间的距离
fit2 <- hclust(d, method="ward.D") #进行Ward层次聚类
plot(fit2) # 绘制树状图展示聚类结果
groups <- cutree(fit2, k=5) # 设定聚类个数为3
# 给聚成的3个类别加上红色边框
rect.hclust(fit2, k=5, border="red")
clustered_data <- data.frame(DataPoint = 1:nrow(loess_predict_mean), Cluster = groups)
write.csv(clustered_data,"/data/CaseControl/S6_cluster_data_AX_9.csv")


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



cluster_data_AX = clustered_data[1:23,]
cluster_data_BIP = clustered_data[24:31,]
cluster_data_DEP = clustered_data[32:81,]
cluster_data_SCH = clustered_data[82:95,]
write.csv(cluster_data_AX, "/data/CaseControl/S6_cluster_data_AX.csv")
write.csv(cluster_data_BIP, "/data/CaseControl/S6_cluster_data_BIP.csv")
write.csv(cluster_data_DEP, "/data/CaseControl/S6_cluster_data_DEP.csv")
write.csv(cluster_data_SCH, "/data/CaseControl/S6_cluster_data_SCH.csv")
#write.csv(cluster_data_SLP, "/Users/yujia/Desktop/a/project/Evolution/MD/MD_prepost/zdata_trajectory/significant_variable/SLP/cluster_data_SLP.csv", row.names = FALSE)
#write.csv(cluster_data_SU, "/Users/yujia/Desktop/a/project/Evolution/MD/MD_prepost/zdata_trajectory/significant_variable/SU/cluster_data_SU.csv", row.names = FALSE)

