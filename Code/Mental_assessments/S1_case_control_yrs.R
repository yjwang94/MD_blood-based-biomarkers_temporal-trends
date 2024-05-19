fpath = "/Users/yujia/Desktop/a/project/Evolution/MD/Data/"
MD = "BIP"
S0_case_control_match_df = read.csv(paste(fpath,"CaseControl/",MD,"/S0_case_control_match_online_followup.csv",sep=''))
S0_TargetData = read.csv("/data/Mental_assessment/S0_TargetData_date_followup.csv")
colnames(S0_case_control_match_df)[1] = "eid"
case = S0_TargetData[S0_TargetData$eid %in% S0_case_control_match_df$eid,]
case_BLFL2MD_yrs = case[c("eid",paste("FL2",MD,"_yrs",sep = ""))]

case_control_BLFL2MD_yrs = merge(S0_case_control_match_df,case_BLFL2MD_yrs,by="eid")

S1_case_control_match_yrs = case_control_BLFL2MD_yrs[,c(1,7)]
S1_case_control_match_yrs$group = 1
for (i in c(2:6)) {
  control = case_control_BLFL2MD_yrs[,c(i,7)]
  colnames(control)[1] = "eid"
  control$group = 0
  S1_case_control_match_yrs = rbind(S1_case_control_match_yrs,control)
}
S1_case_control_match_yrs = S1_case_control_match_yrs[complete.cases(S1_case_control_match_yrs$eid),]
S1_case_control_match_yrs = S1_case_control_match_yrs[complete.cases(S1_case_control_match_yrs$FL2BIP_yrs),]
write.csv(S1_case_control_match_yrs,paste(fpath,'CaseControl/',MD,'/S1_case_control_match_yrs_online_date.csv',sep=""),row.names = FALSE)
