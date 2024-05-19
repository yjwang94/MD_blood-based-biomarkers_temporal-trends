# MD data preparation:DEP,AX,SU,SD

# S0_TargetData to S0_DEP_TargetData
S0_TargetData = read.csv("/data/raw_data/S0_TargetData.csv")

mental_online_date = read.csv("/data/raw_data/online_followup_date.csv")
mental_online_date1 = mental_online_date[c(1,2)]
colnames(mental_online_date1) = c("eid", "mental_online_date")

S0_TargetData_date1 = merge(mental_online_date1, S0_TargetData, by = "eid")

dates <- c(S0_TargetData_date1$AX_date)
S0_TargetData_date1$AX_date = as.Date(dates, "%d/%m/%Y")
S0_TargetData_date1$FL2AX_yrs = (as.Date(S0_TargetData_date1$AX_date)-as.Date(S0_TargetData_date1$mental_online_date))/365
S0_TargetData_date1$FL2AX_yrs = as.numeric(S0_TargetData_date1$FL2AX_yrs)

dates <- c(S0_TargetData_date1$BIP_date)
S0_TargetData_date1$BIP_date = as.Date(dates, "%d/%m/%Y")
S0_TargetData_date1$FL2BIP_yrs = (as.Date(S0_TargetData_date1$BIP_date)-as.Date(S0_TargetData_date1$mental_online_date))/365
S0_TargetData_date1$FL2BIP_yrs = as.numeric(S0_TargetData_date1$FL2BIP_yrs)

dates <- c(S0_TargetData_date1$DEP_date)
S0_TargetData_date1$DEP_date = as.Date(dates, "%d/%m/%Y")
S0_TargetData_date1$FL2DEP_yrs = (as.Date(S0_TargetData_date1$DEP_date)-as.Date(S0_TargetData_date1$mental_online_date))/365
S0_TargetData_date1$FL2DEP_yrs = as.numeric(S0_TargetData_date1$FL2DEP_yrs)

dates <- c(S0_TargetData_date1$SCH_date)
S0_TargetData_date1$SCH_date = as.Date(dates, "%d/%m/%Y")
S0_TargetData_date1$FL2SCH_yrs = (as.Date(S0_TargetData_date1$SCH_date)-as.Date(S0_TargetData_date1$mental_online_date))/365
S0_TargetData_date1$FL2SCH_yrs = as.numeric(S0_TargetData_date1$FL2SCH_yrs)

write.csv(S0_TargetData_date1,"/data/raw_data/S0_TargetData_date_online_followup.csv",row.names = FALSE)
