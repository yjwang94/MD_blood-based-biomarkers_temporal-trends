library("MatchIt")
MD = "BIP"
fpath = "/data/"
S0_TargetData = read.csv("/data/raw_data/S0_TargetData_date_online_followup.csv")
#提取case
cov = read.csv(paste(fpath,"raw_data/Covariates_Data_Imputed.csv",sep = ""))
data_cov = cov
data_cov = merge(data_cov, S0_TargetData[c("eid",paste(MD,"_y",sep=""))], by = "eid")
colnames(data_cov)[13] = "treat"
#rownames(data_cov) = data_cov$eid
case = data_cov[data_cov$treat == 1,]
#提取control
control = data_cov[data_cov$treat == 0,] 
#提取BL2MD_yrs
case_FL2MD_yrs = S0_TargetData[c("eid",paste("FL2",MD,"_yrs",sep = ""))]
case = merge(case,case_FL2MD_yrs,by="eid")
#提取发病前后10年以内的患者
case = case[complete.cases(case[14]),]
case10 = case[case[paste("FL2",MD,"_yrs",sep = "")]<=10 & case[,paste("FL2",MD,"_yrs",sep = "")]>=-10,]
#control组中剔除其他疾病的人群
Blood_immune = read.csv(paste(fpath,'disease_need_exclusion_data/','Blood and immune disorders.csv',sep = ""))
Blood_immune_eid_lst = Blood_immune$eid[Blood_immune$target_y == 1 & Blood_immune$BL2Target_yrs <0]
Mental_behavioural = read.csv(paste(fpath,'disease_need_exclusion_data/','Mental and behavioural disorders.csv',sep = ""))
Mental_behavioural_eid_lst = Mental_behavioural$eid[Mental_behavioural$target_y == 1 & Mental_behavioural$BL2Target_yrs <0]
Infections = read.csv(paste(fpath,'disease_need_exclusion_data/','Diseases of Infections.csv',sep = ""))
Infections_eid_lst = Infections$eid[Infections$target_y == 1 & Infections$BL2Target_yrs <0]
Nervous = read.csv(paste(fpath,'disease_need_exclusion_data/','Nervous system disorders.csv',sep = ""))
Nervous_eid_lst = Nervous$eid[Nervous$target_y == 1 & Nervous$BL2Target_yrs <0]
other_disease = unique(c(Blood_immune_eid_lst,Mental_behavioural_eid_lst,Infections_eid_lst,Nervous_eid_lst ))
control_noohter = control[!control$eid %in% other_disease,]
control_noohter[paste("FL2",MD,"_yrs",sep="")] = 0
all_individuals = rbind(case10,control_noohter)
all_individuals = all_individuals[order(all_individuals$eid),]

m.out1 <- matchit(treat~ + Education + TD_Index + BMI,
                  method = "nearest", distance = "mahalanobis" , ratio = 5,  link = "logit" ,
                  exact = ~ Age + Gender + Ethnicity + SMK_Status + ALC_Status,     #+Smoking+Alcohol+ centre
                  data = all_individuals)

summary(m.out1)
control_row = as.data.frame(m.out1$match.matrix)
control_eid = control_row

for (i in c(1:5)) {
  for (j in c(1:length(control_eid[,i]))) {
    control_eid[j,i] = all_individuals[control_eid[j,i],1]
  }
}

case_control_eid = cbind(case10["eid"],control_eid)

colnames(case_control_eid) = c("case_id","control_id1", "control_id2",	"control_id3","control_id4","control_id5")
for (i in c(1:6)) {
  case_control_eid[,i] = as.character(case_control_eid[,i])
}

#把匹配不到的个体删去
case_control_eid = case_control_eid[complete.cases(case_control_eid$control_id1),]
write.csv(case_control_eid, paste(fpath,'Mental_assessments/',MD,'/S0_case_control_match_online_followup.csv',sep=""),row.names = FALSE)
