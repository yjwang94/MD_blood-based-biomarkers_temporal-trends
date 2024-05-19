install.packages("MatchIt")
library("MatchIt")
#MD = c("Anxiety", "Bipolar", "Depression", "Schizophrenia")
MD = "Depression"
fpath = "/data/"
S0_TargetData = read.csv("/data/raw_data/S0_TargetData.csv")
#case
cov = read.csv(paste(fpath,"Covariates_Data_Imputed.csv",sep = ""))
data_cov = cov
#rownames(data_cov) = data_cov$eid
S0_TargetData1 = S0_TargetData[S0_TargetData[paste(MD,"_date",sep="")] != "",]
data_cov$treat = 0
data_cov$treat[data_cov$eid %in% S0_TargetData1$eid] = 1
case = data_cov[data_cov$treat == 1,]
#control
control = data_cov[data_cov$treat == 0,] 
#BL2MD_yrs
case_BL2MD_yrs = S0_TargetData1[c("eid",paste("BL2",MD,"_yrs",sep = ""))]
case = merge(case,case_BL2MD_yrs,by="eid")

#Patients within 10 years before and after the diagnosis of the disease were extracted
case10 = case[case[paste("BL2",MD,"_yrs",sep = "")]<=10 & case[,paste("BL2",MD,"_yrs",sep = "")]>=-10,]

#Patients with other diseases were excluded from the control group
Blood_immune = read.csv(paste(fpath,'raw_data/disease_need_exclusion_data/','Blood and immune disorders.csv',sep = ""))
Blood_immune_eid_lst = Blood_immune$eid[Blood_immune$target_y == 1 & Blood_immune$BL2Target_yrs <0]
Mental_behavioural = read.csv(paste(fpath,'raw_data/disease_need_exclusion_data/','Mental and behavioural disorders.csv',sep = ""))
Mental_behavioural_eid_lst = Mental_behavioural$eid[Mental_behavioural$target_y == 1 & Mental_behavioural$BL2Target_yrs <0]
Infections = read.csv(paste(fpath,'raw_data/disease_need_exclusion_data/','Diseases of Infections.csv',sep = ""))
Infections_eid_lst = Infections$eid[Infections$target_y == 1 & Infections$BL2Target_yrs <0]
Nervous = read.csv(paste(fpath,'raw_data/disease_need_exclusion_data/','Nervous system disorders.csv',sep = ""))
Nervous_eid_lst = Nervous$eid[Nervous$target_y == 1 & Nervous$BL2Target_yrs <0]
other_disease = unique(c(Blood_immune_eid_lst,Mental_behavioural_eid_lst,Infections_eid_lst,Nervous_eid_lst ))
control_noohter = control[!control$eid %in% other_disease,]
control_noohter[paste("BL2",MD,"_yrs",sep="")] = 0
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

#Individuals that did not match were removed
case_control_eid = case_control_eid[complete.cases(case_control_eid$control_id1),]

write.csv(case_control_eid, paste(fpath,'CaseControl/',MD,'/S0_case_control_match.csv',sep=""),row.names = FALSE)

S0_case_control_match_df = read.csv(paste(fpath,"CaseControl/",MD,"/S0_case_control_match.csv",sep=''))
colnames(S0_case_control_match_df)[1] = "eid"
case = S0_TargetData[S0_TargetData$eid %in% S0_case_control_match_df$eid,]
case_BL2MD_yrs = case[c("eid",paste("BL2",MD,"_yrs",sep = ""))]
case_control_BL2MD_yrs = merge(S1_case_control_match_df,case_BL2MD_yrs,by="eid")

S1_case_control_match_yrs = case_control_BL2MD_yrs[,c(1,7)]
S1_case_control_match_yrs$group = 1
for (i in c(2:6)) {
  control = case_control_BL2MD_yrs[,c(i,7)]
  colnames(control)[1] = "eid"
  control$group = 0
  S1_case_control_match_yrs = rbind(S1_case_control_match_yrs,control)
}
S1_case_control_match_yrs = S1_case_control_match_yrs[complete.cases(S1_case_control_match_yrs$eid),]

write.csv(S1_case_control_match_yrs,paste(fpath,'CaseControl/',MD,'/S1_case_control_match_yrs.csv',sep=""),row.names = FALSE)
