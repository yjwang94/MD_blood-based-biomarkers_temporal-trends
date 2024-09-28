#提取脑区数据
library(stringr)
fpath = "/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/Neuroimaging/"
volume = read.csv(paste(fpath,"Freesurfer_DKT_bl_data.csv",sep = ""))
DTI_field = read.csv(paste(fpath,"Diffusion_MRI_weighted_means.csv",sep = ""))
DTI = read.csv(paste(fpath,"Diffusion_MRI_bl_data.csv",sep = ""))


MD = "SCH"
dpath = "/Users/yujia/Desktop/a/project/Evolution/MD/MD_review/data/CaseControl/"
disease = read.csv(paste(dpath,MD,"/S0_case_control_match.csv",sep = ""))
disease_DKT = volume[volume$eid %in% disease$case_id,]
disease_DKT_DTI = DTI[DTI$eid %in% disease_DKT$eid,]

#合并脑区和cov和显著的blood