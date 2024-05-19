library(dplyr)
library(Hmisc)
#install.packages("Hmisc")
#results从头开始的人口统计学信息统计
#

data_cov = read.csv("/Users/yujia/Desktop/a/project/Evolution/MD/Data/Covariates_Data_Imputed.csv")
head(data_cov)

target = 'AX'
datapath_MD = paste("/Users/yujia/Desktop/a/project/Evolution/MD/Data/CaseControl/",target,"/S1_case_control_match_yrs.csv",sep = '')
data_MD = read.csv(datapath_MD)
data_MD_cov = merge(data_MD, data_cov, by = "eid")


MD = c("BIP", "DEP", "SCH")                                             
target = "AX"
datapath_MD = paste("/Users/yujia/Desktop/a/project/Evolution/MD/Data/CaseControl/",target,"/S1_case_control_match_yrs.csv",sep = '')
MD_all = read.csv(datapath_MD)
#data = read.csv("")
#MD_all = data.frame()
for (target in MD) {
  datapath_MD = paste("/Users/yujia/Desktop/a/project/Evolution/MD/Data/CaseControl/",target,"/S1_case_control_match_yrs.csv",sep = '')
  data = read.csv(datapath_MD)
  MD_all = full_join(MD_all, data, by = "eid")
}
table(is.na(MD_all))

for (i in c(2:length(MD_all))){
  MD_all[is.na(MD_all[,i]),i]=0
}
MD_all$group = 0
for (i in c(1:length(rownames(MD_all)))) {
  if(MD_all$group.x[i] == 1){
    MD_all$group[i] = 1
  }else if(MD_all$group.y[i] == 1){
    MD_all$group[i] = 1
  }else if(MD_all$group.x.x[i] == 1){
    MD_all$group[i] = 1
  }else if(MD_all$group.y.y[i] == 1){
    MD_all$group[i] = 1
  #}else if(MD_all$case_control.x.x.x[i] == 1){
   # MD_all$case_control[i] = 1
  #}else if(MD_all$case_control.y.y.y[i] == 1){
    #MD_all$case_control[i] = 1
  }else{
    MD_all$group[i] = 0
  }
}
MD_all_cov = merge(MD_all, data_cov, by = "eid")


data_MD_cov = MD_all_cov

####### Number #######
number <- data.frame(group=c(), count=c())
for(i in (0:1)){
  subgrp <- data_MD_cov[which(data_MD_cov$group == i), 'group']
  count <- length(subgrp) 
  number <- rbind(number, 
                  data.frame(group = c(i),
                             count = c(count)))
}
number

# Age/TD_index/BMI
summary(data_MD_cov[data_MD_cov$group == 1,]) #case
summary(data_MD_cov[data_MD_cov$group == 0,]) #control

##### Gender #####
#count
sex <- data.frame(group=c(), female=c(), male = c(), ratio = c())
for(i in (0:1)){
  subgrp <- data_MD_cov[which(data_MD_cov$group == i), ]
  sexgrp0 <- subgrp[which(subgrp$Gender == 0),'group']
  sexgrp1 <- subgrp[which(subgrp$Gender == 1),'group']
  female <- length(sexgrp0)
  male <- length(sexgrp1)
  ratio <- male/(female+male)
  sex <- rbind(sex, data.frame(group = c(i),
                               female = c(female),
                               male = c(male),
                               ratio = c(ratio)))
}
sex

#### Ethnicity ####
#count 
Ethnicity = data.frame(group=c(), White = c(), White_ratio = c(), Asian = c(), Asian_ratio = c(), 
                       Black = c(), Black_ratio = c(), Others = c(), Others_ratio = c())
for(i in (0:1)){
  subgrp <- data_MD_cov[which(data_MD_cov$group == i), ]
  grp0 <- subgrp[which(subgrp$Ethnicity == 0),'group']
  grp1 <- subgrp[which(subgrp$Ethnicity == 1),'group']
  grp2 <- subgrp[which(subgrp$Ethnicity == 2),'group']
  grp3 <- subgrp[which(subgrp$Ethnicity == 3),'group']
  White <- length(grp0)
  White_ratio <- White/length(subgrp[,1])
  Black <- length(grp1)
  Black_ratio <- Black/length(subgrp[,1])
  Asian <- length(grp2)
  Asian_ratio <- Asian/length(subgrp[,1])
  Others <- length(grp3)
  Others_ratio <- Others/length(subgrp[,1])
  Ethnicity <- rbind(Ethnicity, data.frame(group = c(i),
                                           White = c(White),
                                           White_ratio = c(White_ratio),
                                           Asian = c(Asian),
                                           Asian_ratio = c(Asian_ratio),
                                           Black = c(Black),
                                           Black_ratio = c(Black_ratio),
                                           Others = c(Others),
                                           Others_ratio = c(Others_ratio)))
}
Ethnicity

### characteristic of mental health assessment group
MD = c("BIP", "DEP", "SCH")                                             
target = "AX"
Neuroticism = read.csv(paste("/Users/yujia/Desktop/a/project/Evolution/MD/MD1205/zdata_trajectory/Neuroticism/zdata_",target,".csv",sep = ""))
Neuroticism_all = Neuroticism
#data = read.csv("")
#MD_all = data.frame()
for (target in MD) {
  Neuroticism = read.csv(paste("/Users/yujia/Desktop/a/project/Evolution/MD/MD1205/zdata_trajectory/Neuroticism/zdata_",target,".csv",sep = ""))
  Neuroticism_all = rbind(Neuroticism_all[,c('eid','group','Neuroticism')], Neuroticism[,c('eid','group','Neuroticism')])
}
table(is.na(Neuroticism_all))
Neuroticism_all_lst = unique(Neuroticism_all$eid[Neuroticism_all$group == 1])

MD = c("BIP", "DEP", "SCH")                                             
target = "AX"
PHQ4 = read.csv(paste("/Users/yujia/Desktop/a/project/Evolution/MD/MD1205/zdata_trajectory/PHQ_4/zdata_",target,".csv",sep = ""))
PHQ4_all = PHQ4
#data = read.csv("")
#MD_all = data.frame()
for (target in MD) {
  PHQ4 = read.csv(paste("/Users/yujia/Desktop/a/project/Evolution/MD/MD1205/zdata_trajectory/PHQ_4/zdata_",target,".csv",sep = ""))
  PHQ4_all = rbind(PHQ4_all[,c('eid','group','PHQ4')], PHQ4[,c('eid','group','PHQ4')])
}
table(is.na(PHQ4_all))
PHQ4_all_lst = unique(PHQ4_all$eid[PHQ4_all$group == 1])


##################### Table1 ####################


target = 'AX'
datapath_MD = paste("/Users/yujia/Desktop/a/project/Evolution/MD/Data/CaseControl/",target,"/S1_case_control_match_yrs.csv",sep = '')
data_MD = read.csv(datapath_MD)
data_cov = read.csv("/Users/yujia/Desktop/a/project/Evolution/MD/Data/Covariates_Data_Imputed.csv")
data_MD_cov = merge(data_MD, data_cov, by = "eid")


####### Number #######
number <- data.frame(group=c(), count=c())
for(i in (0:1)){
  subgrp <- data_MD_cov[which(data_MD_cov$group == i), 'group']
  count <- length(subgrp) 
  number <- rbind(number, 
                  data.frame(group = c(i),
                             count = c(count)))
}
number

# Age/TD_index/BMI
summary(data_MD_cov[data_MD_cov$group == 1,]) #case
summary(data_MD_cov[data_MD_cov$group == 0,]) #control

##### Gender #####
#count
sex <- data.frame(group=c(), female=c(), male = c(), ratio = c())
for(i in (0:1)){
  subgrp <- data_MD_cov[which(data_MD_cov$group == i), ]
  sexgrp0 <- subgrp[which(subgrp$Gender == 0),'group']
  sexgrp1 <- subgrp[which(subgrp$Gender == 1),'group']
  female <- length(sexgrp0)
  male <- length(sexgrp1)
  ratio <- male/(female+male)
  sex <- rbind(sex, data.frame(group = c(i),
                               female = c(female),
                               male = c(male),
                               ratio = c(ratio)))
}
sex

  

#### Ethnicity ####
#count 
Ethnicity = data.frame(group=c(), White = c(), White_ratio = c(), Asian = c(), Asian_ratio = c(), 
                       Black = c(), Black_ratio = c(), Others = c(), Others_ratio = c())
for(i in (0:1)){
  subgrp <- data_MD_cov[which(data_MD_cov$group == i), ]
  grp0 <- subgrp[which(subgrp$Ethnicity == 0),'group']
  grp1 <- subgrp[which(subgrp$Ethnicity == 1),'group']
  grp2 <- subgrp[which(subgrp$Ethnicity == 2),'group']
  grp3 <- subgrp[which(subgrp$Ethnicity == 3),'group']
  White <- length(grp0)
  White_ratio <- White/length(subgrp[,1])
  Black <- length(grp1)
  Black_ratio <- Black/length(subgrp[,1])
  Asian <- length(grp2)
  Asian_ratio <- Asian/length(subgrp[,1])
  Others <- length(grp3)
  Others_ratio <- Others/length(subgrp[,1])
  Ethnicity <- rbind(Ethnicity, data.frame(group = c(i),
                                           White = c(White),
                                           White_ratio = c(White_ratio),
                                           Asian = c(Asian),
                                           Asian_ratio = c(Asian_ratio),
                                           Black = c(Black),
                                           Black_ratio = c(Black_ratio),
                                           Others = c(Others),
                                           Others_ratio = c(Others_ratio)))
}
Ethnicity

#### Education ####
sum(is.na(data_MD_cov$Education))
education = data.frame(group=c(), num0=c(), num0_ratio=c(), num1 = c(), num1_ratio = c(), 
                       num2 = c(), num2_ratio = c())
for(i in (0:1)){
  subgrp <- data_MD_cov[which(data_MD_cov$group == i), ]
  grp0 <- subgrp[which(subgrp$Education <= 2),'group']
  grp1 <- subgrp[which(subgrp$Education <= 6 & subgrp$Education >= 3),'group']
  grp2 <- subgrp[which(subgrp$Education > 6),'group']
  num0 <- length(grp0)
  num0_ratio <- num0/length(subgrp[,1])
  num1 <- length(grp1)
  num1_ratio <- num1/length(subgrp[,1])
  num2 <- length(grp2)
  num2_ratio <- num2/length(subgrp[,1])
  education <- rbind(education, data.frame(group = c(i),
                                           num0 = c(num0),
                                           num0_ratio = c(num0_ratio),
                                           num1 = c(num1),
                                           num1_ratio = c(num1_ratio),
                                           num2 = c(num2),
                                           num2_ratio = c(num2_ratio)))
}
education

#### The numbers of year ####
#count
num_yrs = data.frame(group=c(), num0=c(), num0_ratio=c(), num1 = c(), num1_ratio = c(), 
                     num2 = c(), num2_ratio = c(), num3 = c(), num3_ratio = c())
for(i in (0:1)){
  subgrp <- data_MD_cov[which(data_MD_cov$group == i), ]
  subgrp[,2] = subgrp[,2]*(-1)
  grp0 <- subgrp[which(subgrp[,2] <= -5),'group']
  grp1 <- subgrp[which(subgrp[,2] <= 0 & subgrp[,2] > -5),'group']
  grp2 <- subgrp[which(subgrp[,2] <= 5 & subgrp[,2] > 0),'group']
  grp3 <- subgrp[which(subgrp[,2] > 5),'group']
  num0 <- length(grp0)
  num0_ratio <- num0/length(subgrp[,1])
  num1 <- length(grp1)
  num1_ratio <- num1/length(subgrp[,1])
  num2 <- length(grp2)
  num2_ratio <- num2/length(subgrp[,1])
  num3 <- length(grp3)
  num3_ratio <- num3/length(subgrp[,1])
  num_yrs <- rbind(num_yrs, data.frame(group = c(i),
                                           num0 = c(num0),
                                           num0_ratio = c(num0_ratio),
                                           num1 = c(num1),
                                           num1_ratio = c(num1_ratio),
                                           num2 = c(num2),
                                           num2_ratio = c(num2_ratio),
                                           num3 = c(num3),
                                           num3_ratio = c(num3_ratio)))
}
num_yrs


