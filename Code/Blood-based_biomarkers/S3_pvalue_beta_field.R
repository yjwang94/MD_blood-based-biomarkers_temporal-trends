
UKB_pheno = read.csv("/data/PhenoData2Analysis/UKB_pheno.csv")
UKB_blood_dict = read.csv("/data/PhenoData2Analysis/UKB_Blood_dict.csv")
NMR_dict = read.csv("/data/PhenoData2Analysis/abbreviation_NMR.csv")
Serum_dict = read.csv("/data/PhenoData2Analysis/Serum_dict1.csv")

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

#将NMR_dict加上FieldID
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
colnames(data_FieldID)[1] = "FieldID_full"

MD = "Depression"
fpath = "/data/CaseControl/"
data = read.csv(paste(fpath,MD,'/S3_Pvalue_beta.csv',sep = ""))
data = data[,c("FieldID_full","beta_group","p_group_bfi","beta_group_time","p_group_time_bfi","beta_group_time2","p_group_time2_bfi")]
data_FieldID_beta_p = merge(data_FieldID,data, by = 'FieldID_full')
colnames(UKB_blood_dict)[1] = "FieldID_full"

write.csv(data_FieldID_beta_p,paste(fpath,MD,'/S3_Pvalue_beta_field.csv',sep = ""),row.names = FALSE)

