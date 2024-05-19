## Main workflow and corresponding Rscript

### The temporal trend of blood-based biomarkers

#### 1. Match the corresponding controls for the case group (S1_case_control_matched.R)

For each participant with a mental disorder, nearly 5 controls were selected based on propensity score matching for age, sex, ethnicity, centers, years of education, body mass index (BMI), socioeconomical status, drinking status and smoking status using the R package **“Matchit”**. Controls were defined as individuals without blood and immune disorders, mental and behavioral disorders, disease of infections and nervous system disorders. We used the baseline date and the diagnosis date to calculate the years from clinical diagnosis (date of clinical diagnosis minus baseline date). The time since diagnosis of control participants were arranged as the time since diagnosis of their matched patients. Consequently, each patient and matched control participants had the same index date. The precessed data was recorded in "S0_case_control_match_yrs.csv". Post-match data for all four diseases were stored in the folder "CaseControl".



#### 2.Prepare phenotypic data (S2_pheno_prepare.R)

In this section, we prepare 227 blood-based biomarkers for analysis. 

#### 3.Establish multiple linear regression models(S3_multiple_regression & S3_pvalue_beta_field.R)

Our primary objective was to identify measures of divergence evolution between individuals with mental disorders and controls. We used a multivariable linear regression model to compare each variable of interest between the two groups. We included group and its interactions with time and the square of time as independent variables, with blood-based markers as the response variable. Other independent variables included age, sex, ethnicity, years of education, townsend deprivation index, study site, BMI, drinking status, and smoking status at the baseline visit. Ethnicity was categorized as white, asian, black, and other ethnicity (e.g., mixed ethnicity or other ethnicity). Continuous variables were normalized before regression, and categorical variables were made dummy variables. Missing values were not included in the analysis. P-values for group and its interactions with time and the square of time were corrected for multiple comparisons using the Bonferroni method, respectively. Blood-based markers with any corrected p-values less than 0.01 were used for temporal trends visualization. The analysis was performed using the **StatsModel** (v0.11.1) and **ScikitLearn** (v0.24.1) packages in **Python** (v3.9).The result of multiple linear regression model was recorded in file "S3_Pvalue_beta.csv".

#### 4. Temporal trends visualization(S4_zdata_prepare.R, S5_loess_predict_mean.R, S6_clusterdata_prepare.R, S7_zdata_trajectory.R)

To visualize the temporal trends of blood-based markers in mental disorders, we employed the locally estimated scatterplot smoothing (LOESS) method to assess linear or non-linear patterns in blood-based marker evolution over time. To ensure that each assessment's values were easily comparable, blood-based marker values for mental disorders were transformed using z-scores relative to controls, adjusting for the mentioned covariates. To reduce the complexity of blood-based markers, we grouped blood-based markers with similar temporal trends across diseases using unsupervised hierarchical clustering via the hclust function from the R stats package. 

### The temporal trend of  mental assessments

The analysis of the temporal trend of  mental assessments was recorded in file "/Code/Mental_assessments/".

### Correlation between the identified blood-based markers and brain structures

The correlation between the identified blood-based markers and brain structures were performed using Pearson's correlation adjusted for the aforementioned covariates.



## File description

| File name | Description                           |
| --------- | ------------------------------------- |
| code      | data process and result analysis code |
| Figure    | The Figure of the analysis result     |
| ReadMe.md | ReadMe Markdown                       |





