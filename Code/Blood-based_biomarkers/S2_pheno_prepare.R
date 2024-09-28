library(ggplot2)
library(reshape2)
library(dplyr)

#MD = c("Anxiety", "Bipolar", "Depression", "Schizophrenia")
MD = "Depression"
fpath = paste('/data/CaseControl/', MD, '/S1_case_control_match_yrs.csv', sep = '')
data = read.csv(fpath)
covariance = read.csv("/data/raw_data//Covariates_Data_Imputed.csv")
UKB_pheno = read.csv("/data/PhenoData2Analysis/UKB_pheno.csv")
UKB_blood_dict = read.csv("/data/PhenoData2Analysis/UKB_Blood_dict.csv")
NMR_dict = read.csv("/U/data/PhenoData2Analysis/abbreviation_NMR.csv")
Serum_dict = read.csv("/data/PhenoData2Analysis/Serum_dict1.csv")

#zdata
data_cov = merge(data, covariance, by = 'eid')

# rename pheno of UKB_blood_dict
UKB_blood_dict$FieldID = paste(UKB_blood_dict$FieldID, '-0.0',sep = '')
Serum_dict$FieldID = paste(Serum_dict$FieldID,'-0.0', sep = '')
UKB_blood_dict$Field = tolower(UKB_blood_dict$Field)
NMR_dict$Description = tolower(NMR_dict$Description)
Serum_dict$Field = tolower(Serum_dict$Field)

# rename pheno of UKB_pheno
tmp = colnames(UKB_pheno)[2:383]
for (i in (1:382)) {
  tmp[i] = sub('X', '', tmp[i])
  tmp[i] = sub('\\.', '\\-', tmp[i])
}
colnames(UKB_pheno)[2:383] = tmp

#Add FieldID
for (i in NMR_dict$Description) {
  j = UKB_blood_dict$FieldID[UKB_blood_dict$Field == i]
  NMR_dict$FieldID[NMR_dict$Description == i] = j
}

data_NMR = UKB_blood_dict[UKB_blood_dict$Main_Category == "Metabolomics",]
for (i in data_NMR$Field) {
  j = NMR_dict$Abbreviation[NMR_dict$Description == i]
  UKB_blood_dict$field[UKB_blood_dict$Field == i] = j
}

#Serum_dict$Sub_Category1 = tolower(Serum_dict$Sub_Category1)

data_Serum = UKB_blood_dict[UKB_blood_dict$Main_Category == "Serum",]
for (i in data_Serum$FieldID) {
  j = Serum_dict$Sub_Category1[Serum_dict$FieldID == i]
  UKB_blood_dict$Sub_Category[UKB_blood_dict$FieldID == i] = j
}
for (i in data_Serum$FieldID) {
  j = Serum_dict$Abbreviation[Serum_dict$FieldID == i]
  UKB_blood_dict$field[UKB_blood_dict$FieldID == i] = j
}


Serum_Field = Serum_dict[,c(1, 6)]
NMR_Field = NMR_dict[,c(8,1)]
data_FieldID = rbind(Serum_Field, NMR_Field)


#UKB_pheno227
UKB_pheno227 = UKB_pheno[, c('eid', UKB_blood_dict$FieldID)]

data_cov_pheno = merge(data_cov, UKB_pheno227, by = 'eid')
colnames(data_cov_pheno)[2] = "duration"
outfile = paste('/data/CaseControl/',MD,'/S2_data_cov_pheno_forPvalue.csv',sep = '')
write.csv(data_cov_pheno, outfile,row.names = FALSE)
