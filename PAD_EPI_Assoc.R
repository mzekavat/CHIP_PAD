### PAD Project

library(data.table)
library(survival)
library(tableone)

#########################################################################################################
### UKBB Phenos & QC:

pheno = fread("/medpop/esp2/mzekavat/UKBB/ukbb_PhenoFile.ALL_500k.UpdatedIncdPhenos_202020.txt")
pheno = data.frame(pheno)
pheno$IDs_toRemove_SampleQCv2 = ifelse((!(pheno$Submitted_Gender == pheno$Inferred_Gender) |  pheno$Second_deg_relOrHigher_toRemove ==1 | pheno$Non_Consented== 1),1,0) ## N:
pheno = pheno[which(is.na(pheno$IDs_toRemove_SampleQCv2 ) & pheno$inWES == 1),]
pheno = pheno[-which(pheno$Prev_Maryam_HemeCa_Phenos ==1),]
pheno = pheno[-which(pheno$Prev_Peripheral_vascular_disease ==1),]

### ------Baseline Summary Stats:
myvars = c("age", "Sex", "Race","SmokingStatusv2","Alcohol_intake_freq_last4wkNum","Freq_Excercise_last4wkNum","Townsend", "StressInPast2yr", "N_Sweets_Handfulls", "N_Veggi_Servings","BMI", "Prev_Diabetes_Type_2", "Prev_Coronary_Artery_Disease_SOFT", "Prev_Hypertension", "Prev_Hypercholesterolemia")

catVars <- c("Sex", "Race","StressInPast2yr","SmokingStatusv2","Prev_Diabetes_Type_2", "Prev_Coronary_Artery_Disease_SOFT", "Prev_Hypertension", "Prev_Hypercholesterolemia")

tab2 <- CreateTableOne(vars = myvars, data = pheno, factorVars = catVars, strata = c("NewCHIP"))
tab2Mat <- print(tab2, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(tab2Mat, file = "/medpop/esp2/mzekavat/UKBB/PAD_CHIP/UKBBsummaryTable_byCHIP.csv")


myvars = c("age", "Sex", "Race","SmokingStatusv2","BMI", "Prev_Diabetes_Type_2", "Prev_Coronary_Artery_Disease_SOFT", "Prev_Hypertension", "Prev_Hypercholesterolemia")

catVars <- c("Sex", "Race","SmokingStatusv2","Prev_Diabetes_Type_2", "Prev_Coronary_Artery_Disease_SOFT", "Prev_Hypertension", "Prev_Hypercholesterolemia")

tab2 <- CreateTableOne(vars = myvars, data = pheno, factorVars = catVars, strata = c("NewCHIP"))
tab2Mat <- print(tab2, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(tab2Mat, file = "/medpop/esp2/mzekavat/UKBB/PAD_CHIP/UKBBsummaryTable_byCHIP.simpler.csv")

### ------Multivariate association of CHIP with possible risk factors:
pheno$SmokingStatusv2 = factor(pheno$SmokingStatusv2 , levels=c("Never",  "Previous","Current"), ordered=FALSE)
res = summary(glm(NewCHIP ~ Sex + age + age2+ SmokingStatusv2+Alcohol_intake_freq_last4wkNum + Freq_Excercise_last4wkNum + scale(Townsend) + N_Sweets_Handfulls+ scale(BMI) + Prev_Diabetes_Type_2 + Prev_Coronary_Artery_Disease_SOFT + Prev_Hypertension + Prev_Hypercholesterolemia + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data=pheno, family='binomial'))
write.csv(data.frame(res$coeff), file = "/medpop/esp2/mzekavat/UKBB/PAD_CHIP/UKBBriskFactorMultivarAssoc.csv")

pheno$CHIP_DNMT3A = ifelse(pheno$NewCHIP == 0, 0, 
									ifelse(pheno$NewCHIP ==1 & pheno$NewCHIP_Hugo_Symbol %in% c("DNMT3A"), 1,
										ifelse(pheno$NewCHIP ==1 & !(pheno$NewCHIP_Hugo_Symbol %in% c("DNMT3A")),NA, 0)))

pheno$CHIP_TET2 = ifelse(pheno$NewCHIP == 0, 0, 
									ifelse(pheno$NewCHIP ==1 & pheno$NewCHIP_Hugo_Symbol %in% c("TET2"), 1,
										ifelse(pheno$NewCHIP ==1 & !(pheno$NewCHIP_Hugo_Symbol %in% c("TET2")),NA, 0)))

pheno$CHIP_JAK2 = ifelse(pheno$NewCHIP == 0, 0, 
									ifelse(pheno$NewCHIP ==1 & pheno$NewCHIP_Hugo_Symbol %in% c("JAK2") | pheno$V617clonal == 1, 1,
										ifelse(pheno$NewCHIP ==1 & !(pheno$NewCHIP_Hugo_Symbol %in% c("JAK2")) ,NA, 0)))

pheno$CHIP_ASXL1 = ifelse(pheno$NewCHIP == 0, 0, 
									ifelse(pheno$NewCHIP ==1 & pheno$NewCHIP_Hugo_Symbol %in% c("ASXL1"), 1,
										ifelse(pheno$NewCHIP ==1 & !(pheno$NewCHIP_Hugo_Symbol %in% c("ASXL1")),NA, 0)))

pheno$LargeCHIP_DNMT3A = ifelse(pheno$New_Large_CHIP_01 == 0, 0, 
									ifelse(pheno$New_Large_CHIP_01 ==1 & pheno$NewCHIP_Hugo_Symbol %in% c("DNMT3A"), 1,
										ifelse(pheno$New_Large_CHIP_01 ==1 & !(pheno$NewCHIP_Hugo_Symbol %in% c("DNMT3A")),NA, 0)))

pheno$LargeCHIP_TET2 = ifelse(pheno$New_Large_CHIP_01 == 0, 0, 
									ifelse(pheno$New_Large_CHIP_01 ==1 & pheno$NewCHIP_Hugo_Symbol %in% c("TET2"), 1,
										ifelse(pheno$New_Large_CHIP_01 ==1 & !(pheno$NewCHIP_Hugo_Symbol %in% c("TET2")),NA, 0)))

pheno$LargeCHIP_JAK2 = ifelse(pheno$New_Large_CHIP_01 == 0, 0, 
									ifelse(pheno$New_Large_CHIP_01 ==1 & pheno$NewCHIP_Hugo_Symbol %in% c("JAK2") | pheno$V617clonal == 1, 1,
										ifelse(pheno$New_Large_CHIP_01 ==1 & !(pheno$NewCHIP_Hugo_Symbol %in% c("JAK2")) ,NA, 0)))


pheno$LargeCHIP_ASXL1 = ifelse(pheno$New_Large_CHIP_01 == 0, 0, 
									ifelse(pheno$New_Large_CHIP_01 ==1 & pheno$NewCHIP_Hugo_Symbol %in% c("ASXL1"), 1,
										ifelse(pheno$New_Large_CHIP_01 ==1 & !(pheno$NewCHIP_Hugo_Symbol %in% c("ASXL1")),NA, 0)))

pheno$CHIP_minusJAK2 = ifelse(pheno$NewCHIP == 0, 0, 
									ifelse(pheno$NewCHIP ==1 & !pheno$NewCHIP_Hugo_Symbol %in% c("JAK2"), 1,
										ifelse(pheno$NewCHIP ==1 & (pheno$NewCHIP_Hugo_Symbol %in% c("JAK2")),NA, 0)))


pheno$LargeCHIP_minusJAK2 = ifelse(pheno$New_Large_CHIP == "NO_CHIP", "NO_CHIP",
 								ifelse(pheno$New_Large_CHIP == "SMALL_CHIP" & pheno$NewCHIP_Hugo_Symbol %in% c("JAK2"),NA,
									ifelse(pheno$New_Large_CHIP == "SMALL_CHIP" & !(pheno$NewCHIP_Hugo_Symbol %in% c("JAK2")), 'SMALL_CHIP',
										ifelse(pheno$New_Large_CHIP == "LARGE_CHIP" & pheno$NewCHIP_Hugo_Symbol %in% c("JAK2"),NA,
									ifelse(pheno$New_Large_CHIP == "LARGE_CHIP" & !(pheno$NewCHIP_Hugo_Symbol %in% c("JAK2")), 'LARGE_CHIP', NA)))))

pheno$New_Large_CHIP = factor(pheno$New_Large_CHIP , levels=c("NO_CHIP",  "LARGE_CHIP","SMALL_CHIP"), ordered=FALSE)
pheno$LargeCHIP_minusJAK2 = factor(pheno$LargeCHIP_minusJAK2 , levels=c("NO_CHIP",  "LARGE_CHIP","SMALL_CHIP"), ordered=FALSE)


CHIP_Names = c('NewCHIP','New_Large_CHIP')
pheno_list2 = c('Peripheral_vascular_disease')
vars = c(CHIP_Names)
summaryDF = data.frame()
for (i in 1:length(pheno_list2)){
	print(pheno_list2[i])
for (j in 1:length(vars)){
	tmp = pheno
	print(vars[j])
	#if(vars[j] %in%CHIP_Names){tmp = pheno}
	#if(vars[j] %in%c(CHUD_Names, CHUD_Normalized_Names)){tmp = pheno[-which(pheno$hasCHIP ==1),]}
indx = which(colnames(tmp) == vars[j])[1]

#removing prevalent cases
preval_col = which(colnames(tmp) == paste("Prev_",pheno_list2[i],sep=""))
incid_col = which(colnames(tmp) == paste("Incd_",pheno_list2[i],sep=""))

if (length(which(tmp[,preval_col] == 1))>0 ){tmp = tmp[-which(tmp[,preval_col] == 1),]}

if (length(which(tmp[,incid_col] ==1  &(tmp[,indx]>0 | tmp[,indx] == "LARGE_CHIP") & !is.na(tmp$age) & !is.na(tmp$PC1)  & !is.na(tmp$smok_detailed_))) > 4){

#running survival anal
tmp$SurvObj <- Surv(tmp[,which(colnames(tmp) == paste("FollowUp_",pheno_list2[i], sep=""))],  tmp[,which(colnames(tmp) == paste("Incd_",pheno_list2[i], sep=""))]== 1)


df = as.data.frame(t(as.data.frame(summary( coxph(tmp$SurvObj ~ tmp[,indx], data = tmp))$coeff[1,])))
df$Adjustment = "Unadjusted"
df$y = pheno_list2[i]
df$x = vars[j]
df$N_Incd_cases = length(which(tmp[,incid_col] == 1 ))
df$N_Controls = length(which(tmp[,incid_col] == 0 ))
df$N_Incd_cases_withVar = length(which(tmp[,incid_col] == 1 & (tmp[,indx]>0 | tmp[,indx] == "LARGE_CHIP") ))
df$N_Controls_withVar = length(which(tmp[,incid_col] == 0 & (tmp[,indx]>0 | tmp[,indx] == "LARGE_CHIP") ))
summaryDF = rbind(summaryDF, df)

df = as.data.frame(t(as.data.frame(summary( coxph(tmp$SurvObj ~ tmp[,indx]+age +age2 + Sex_numeric+SmokingStatusv2 +scale(Townsend)+ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = tmp))$coeff[1,])))
df$Adjustment = "Sparsely Adjusted"
df$y = pheno_list2[i]
df$x = vars[j]
df$N_Incd_cases = length(which(tmp[,incid_col] == 1  & !is.na(tmp$age) & !is.na(tmp$PC1)  & !is.na(tmp$SmokingStatusv2)& !is.na(tmp$Townsend)))
df$N_Controls = length(which(tmp[,incid_col] == 0 & !is.na(tmp$age) & !is.na(tmp$PC1)  & !is.na(tmp$SmokingStatusv2)& !is.na(tmp$Townsend)))
df$N_Incd_cases_withVar = length(which(tmp[,incid_col] == 1 & (tmp[,indx]>0 | tmp[,indx] == "LARGE_CHIP") &!is.na(tmp$age) & !is.na(tmp$PC1)  & !is.na(tmp$SmokingStatusv2)& !is.na(tmp$Townsend)))
df$N_Controls_withVar = length(which(tmp[,incid_col] == 0 & (tmp[,indx]>0 | tmp[,indx] == "LARGE_CHIP") &!is.na(tmp$age) & !is.na(tmp$PC1)  & !is.na(tmp$SmokingStatusv2)& !is.na(tmp$Townsend)))
summaryDF = rbind(summaryDF, df)


df = as.data.frame(t(as.data.frame(summary( coxph(tmp$SurvObj ~ tmp[,indx]+age +age2 + Sex_numeric+SmokingStatusv2 + scale(Townsend)+Prev_Diabetes_Type_2+Prev_Hypercholesterolemia+Prev_Hypertension+ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = tmp))$coeff[1,])))#Prev_Coronary_Artery_Disease_SOFT
df$Adjustment = "Fully Adjusted"
df$y = pheno_list2[i]
df$x = vars[j]
df$N_Incd_cases = length(which(tmp[,incid_col] == 1  & !is.na(tmp$age) & !is.na(tmp$PC1)  & !is.na(tmp$SmokingStatusv2)& !is.na(tmp$Townsend)))
df$N_Controls = length(which(tmp[,incid_col] == 0 & !is.na(tmp$age) & !is.na(tmp$PC1)  & !is.na(tmp$SmokingStatusv2)& !is.na(tmp$Townsend)))
df$N_Incd_cases_withVar = length(which(tmp[,incid_col] == 1 & (tmp[,indx]>0 | tmp[,indx] == "LARGE_CHIP") &!is.na(tmp$age) & !is.na(tmp$PC1)  & !is.na(tmp$SmokingStatusv2)& !is.na(tmp$Townsend)))
df$N_Controls_withVar = length(which(tmp[,incid_col] == 0 & (tmp[,indx]>0 | tmp[,indx] == "LARGE_CHIP") &!is.na(tmp$age) & !is.na(tmp$PC1)  & !is.na(tmp$SmokingStatusv2)& !is.na(tmp$Townsend)))

summaryDF = rbind(summaryDF, df)
}
}
}

write.table(summaryDF,"/medpop/esp2/mzekavat/UKBB/PAD_CHIP/UKB_PAD_CHIP_assoc.txt",
col.names = T, row.names = F, quote = F, sep = "\t")

### ------ Gene-specific associations:

CHIP_Names = c('CHIP_minusJAK2','LargeCHIP_minusJAK2','CHIP_DNMT3A', 'CHIP_TET2', 'CHIP_JAK2', 'CHIP_ASXL1', 'LargeCHIP_DNMT3A', 'LargeCHIP_TET2', 'LargeCHIP_JAK2', 'LargeCHIP_ASXL1')
pheno_list2 = c('Peripheral_vascular_disease')
vars = c(CHIP_Names)
summaryDF = data.frame()
for (i in 1:length(pheno_list2)){
	print(pheno_list2[i])
for (j in 1:length(vars)){
	tmp = pheno
	print(vars[j])
	#if(vars[j] %in%CHIP_Names){tmp = pheno}
	#if(vars[j] %in%c(CHUD_Names, CHUD_Normalized_Names)){tmp = pheno[-which(pheno$hasCHIP ==1),]}
indx = which(colnames(tmp) == vars[j])[1]

#removing prevalent cases
preval_col = which(colnames(tmp) == paste("Prev_",pheno_list2[i],sep=""))
incid_col = which(colnames(tmp) == paste("Incd_",pheno_list2[i],sep=""))

if (length(which(tmp[,preval_col] == 1))>0 ){tmp = tmp[-which(tmp[,preval_col] == 1),]}

if (length(which(tmp[,incid_col] ==1  &(tmp[,indx]>0 | tmp[,indx] == "LARGE_CHIP") & !is.na(tmp$age) & !is.na(tmp$PC1)  & !is.na(tmp$smok_detailed_))) > 4){

#running survival anal
tmp$SurvObj <- Surv(tmp[,which(colnames(tmp) == paste("FollowUp_",pheno_list2[i], sep=""))],  tmp[,which(colnames(tmp) == paste("Incd_",pheno_list2[i], sep=""))]== 1)


df = as.data.frame(t(as.data.frame(summary( coxph(tmp$SurvObj ~ tmp[,indx]+age +age2 + Sex_numeric+SmokingStatusv2 + scale(Townsend)+Prev_Diabetes_Type_2+Prev_Hypercholesterolemia+Prev_Hypertension+ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = tmp))$coeff[1,])))
df$Adjustment = "Sparsely Adjusted"
df$y = pheno_list2[i]
df$x = vars[j]
df$N_Incd_cases = length(which(tmp[,incid_col] == 1  & !is.na(tmp$age) & !is.na(tmp$PC1)  & !is.na(tmp$SmokingStatusv2)& !is.na(tmp$Townsend)))
df$N_Controls = length(which(tmp[,incid_col] == 0 & !is.na(tmp$age) & !is.na(tmp$PC1)  & !is.na(tmp$SmokingStatusv2)& !is.na(tmp$Townsend)))
df$N_Incd_cases_withVar = length(which(tmp[,incid_col] == 1 & (tmp[,indx]>0 | tmp[,indx] == "LARGE_CHIP") &!is.na(tmp$age) & !is.na(tmp$PC1)  & !is.na(tmp$SmokingStatusv2)& !is.na(tmp$Townsend)))
df$N_Controls_withVar = length(which(tmp[,incid_col] == 0 & (tmp[,indx]>0 | tmp[,indx] == "LARGE_CHIP") &!is.na(tmp$age) & !is.na(tmp$PC1)  & !is.na(tmp$SmokingStatusv2)& !is.na(tmp$Townsend)))
summaryDF = rbind(summaryDF, df)

}
}
}

write.table(summaryDF,"/medpop/esp2/mzekavat/UKBB/PAD_CHIP/UKB_PAD_CHIP_assoc.GeneSpecific.txt",
col.names = T, row.names = F, quote = F, sep = "\t")

##################################################################################################################################################################################################################
### MGBB Phenos:
library(data.table)
library(survival)
library(tableone)
pheno = fread("/medpop/esp2/mzekavat/Partners_BB/phenos/all_basic_MGBB_phenos.plus_mCAs_CHIP.txt")
pheno = data.frame(pheno)
pheno = pheno[-which((pheno$sexCheck_toRemove==1 & !is.na(pheno$sexCheck_toRemove))| (pheno$first_or_SecondDegRelative_toRemove == 1 & !is.na(pheno$first_or_SecondDegRelative_toRemove)) | (pheno$inWES !=1 | is.na(pheno$inWES) ) ),]
pheno = pheno[-which(pheno$Prev_Heme_CA ==1 | pheno$age<20),]
pheno = pheno[-which(pheno$Gender == "U" | pheno$Prev_pad.x== 1),]
pheno$SmokingStatus = ifelse(pheno$Currently_Smoke != "" , "Current Smoker", 
						ifelse(pheno$Smoked_Previously != "" | pheno$Smoking_Final ==1, "Prior Smoker", 
							ifelse(pheno$Never_Smoked != "" | pheno$Smoking_Final == 0, "Never Smoker", "Unknown")))
pheno$Nearest_BMI_to_Collect.Date = ifelse(pheno$Nearest_BMI_to_Collect.Date<10 | pheno$Nearest_BMI_to_Collect.Date>60, NA, pheno$Nearest_BMI_to_Collect.Date)
pheno$age2 = pheno$age*pheno$age

pheno$Large_CHIP = factor(pheno$Large_CHIP , levels=c("NO_CHIP",  "LARGE_CHIP","SMALL_CHIP"), ordered=FALSE)

vars = c( 'age', 'Gender','Race', 'SmokingStatus','Nearest_BMI_to_Collect.Date', 'Prev_Diabetes._Type_2', 'Prev_cad.x', 'Prev_htn.x', 'Prev_hyperlipidemia.x' )


catVars <- c("Gender", "SmokingStatus", "Race", "Prev_Diabetes._Type_2", "Prev_cad.x", "Prev_htn.x", "Prev_hyperlipidemia.x")

tab2 <- CreateTableOne(vars = vars, data = pheno, factorVars = catVars,strata = c("CHIP") )
tab2Mat <- print(tab2, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(tab2Mat, file = "/medpop/esp2/mzekavat/UKBB/PAD_CHIP/MGBBsummaryTable_byCHIP.csv")

pheno = pheno[,-c(262:269)]
colnames(pheno)[198] = "Prev_pad"


pheno_list2 = c('pad')
CHIP_Names = c('CHIP','Large_CHIP')
vars = c(CHIP_Names)
summaryDF = data.frame()
for (i in 1:length(pheno_list2)){
  print(pheno_list2[i])
  tmp = pheno

#removing prevalent cases
  preval_col = which(colnames(tmp) == paste("Prev_",pheno_list2[i],sep=""))
  incid_col = which(colnames(tmp) == paste("Incd_",pheno_list2[i],sep=""))

for (j in 1:length(vars)){
print(vars[j])
indx = which(colnames(tmp) == vars[j])[1]

if (length(which(tmp[,incid_col] ==1  &(tmp[,indx]>0 | tmp[,indx] == "LARGE_CHIP") & !is.na(tmp$age) & !is.na(tmp$Race)  & !is.na(tmp$SmokingStatus))) > 1){

#running survival anal
tmp$SurvObj <- Surv(tmp[,which(colnames(tmp) == paste("FollowUp_",pheno_list2[i],sep=""))],  tmp[,which(colnames(tmp) == paste("Incd_",pheno_list2[i], sep=""))]== 1)

df = as.data.frame(t(as.data.frame(summary( coxph(tmp$SurvObj ~ tmp[,indx]+age +age2 + Gender+SmokingStatus +PC1_WES + PC2_WES + PC3_WES + PC4_WES + PC5_WES + PC6_WES + PC7_WES + PC8_WES + PC9_WES + PC10_WES, data = tmp))$coeff[1,])))
df$y = pheno_list2[i]
df$x = vars[j]
df$N_Incd_cases = length(which(tmp[,incid_col] == 1  & !is.na(tmp$age) & !is.na(tmp$PC1_WES)  & !is.na(tmp$SmokingStatus)))
df$N_Controls = length(which(tmp[,incid_col] == 0 & !is.na(tmp$age) & !is.na(tmp$PC1_WES)  & !is.na(tmp$SmokingStatus)))
df$N_Incd_cases_withVar = length(which(tmp[,incid_col] == 1 & (tmp[,indx]>0 | tmp[,indx] == "LARGE_CHIP") &!is.na(tmp$age) & !is.na(tmp$PC1_WES)  & !is.na(tmp$SmokingStatus)))
df$N_Controls_withVar = length(which(tmp[,incid_col] == 0 & (tmp[,indx]>0 | tmp[,indx] == "LARGE_CHIP") &!is.na(tmp$age) & !is.na(tmp$PC1_WES)  & !is.na(tmp$SmokingStatus)))

summaryDF = rbind(summaryDF, df)
}
}
}

write.table(summaryDF,paste("/medpop/esp2/mzekavat/UKBB/PAD_CHIP/MGBB_CHIP_PADassoc.txt", sep=""),
col.names = T, row.names = F, quote = F, sep = "\t")




#-------- Meta -analyzing PAD association between MGBB and UKBB and making Forest plot:

MGB_Forest = read.table("/medpop/esp2/mzekavat/UKBB/PAD_CHIP/MGBB_CHIP_PADassoc.txt",
                   header=T, as.is=T, stringsAsFactors=F, comment.char = '', sep="\t")
MGB_Forest$Cohort = "MGB Biobank"
UKB_Forest = read.table("/medpop/esp2/mzekavat/UKBB/PAD_CHIP/UKB_PAD_CHIP_assoc.txt",
                   header=T, as.is=T, stringsAsFactors=F, comment.char = '', sep="\t")
UKB_Forest = UKB_Forest[which(UKB_Forest$Adjustment == "Fully Adjusted"),]
UKB_Forest=UKB_Forest[,-c(6)]
UKB_Forest$Cohort = "UK Biobank"

#Forest$Pval = as.character(Forest$Pr...z..)
Forest = rbind(UKB_Forest, MGB_Forest)
Forest$Pval = Forest$Pr...z..


Forest$HR = Forest$exp.coef.
Forest$N_Controls = formatC(Forest$N_Controls, format = "f", big.mark = ",", drop0trailing = TRUE)
Forest$N_Incd_cases = formatC(Forest$N_Incd_cases, format = "f", big.mark = ",", drop0trailing = TRUE)
Forest$N_Incd_cases_withVar = formatC(Forest$N_Incd_cases_withVar, format = "f", big.mark = ",", drop0trailing = TRUE)
Forest$N_Controls_withVar = formatC(Forest$N_Controls_withVar, format = "f", big.mark = ",", drop0trailing = TRUE)
Forest <- Forest[order(Forest$Pr...z.. ),]


Forest$x = ifelse(Forest$x %in% c("CHIP", "NewCHIP"), "CHIP","Large CHIP")
CHIPonly = Forest

library(meta)
library(grid)
library(scales)

CHIPonly <- CHIPonly[order(CHIPonly$Pr...z..,CHIPonly$x  ),]
summarymeta=metagen(coef,se.coef.,studlab=Cohort,byvar = x,data=CHIPonly,sm="HR",comb.fixed = TRUE)

metaAnal = data.frame(Cohort = "Overall", x = summarymeta$bylevs,y = "Peripheral_vascular_disease", coef = summarymeta$TE.fixed.w , se.coef. = summarymeta$seTE.fixed.w, Pval = summarymeta$pval.fixed.w,Pr...z.. = summarymeta$pval.fixed.w, HR = exp(summarymeta$TE.fixed.w))
metaAnal$N_Incd_cases = ""
metaAnal$N_Controls =""
metaAnal$N_Incd_cases_withVar =""
metaAnal$N_Controls_withVar =""
CHIPonly = CHIPonly[,-c(2,4)]
CHIPonly$Pval = CHIPonly$Pr...z..

mergedv2 = rbind(CHIPonly, metaAnal)
mergedv2$Pval = formatC(mergedv2$Pval, format = "g", digits = 2)
mergedv2$x = ordered(mergedv2$x, levels=c("CHIP", "Large CHIP"))
mergedv2$Cohort = ordered(mergedv2$Cohort, levels=c('UK Biobank', 'MGB Biobank', 'Overall'))
mergedv2 = mergedv2[order(mergedv2$x, mergedv2$Cohort),]

summarymeta=metagen(coef,se.coef.,studlab=Cohort,byvar = x,data=mergedv2,sm="HR")

pdf(paste("/medpop/esp2/mzekavat/UKBB/PAD_CHIP/Forest_CHIP.PAD.pdf",sep=""), width = 13, height= 3.5)
forest(summarymeta,leftcols =c("studlab"),leftlabs = c("          "),just='left',colgap.left=unit(0.1,"cm"),xlim=c(0.8, 6),at=c(0.8,1,2,4,6),
      smlab='Association of CHIP\nwith Incident Peripheral Artery Disease',colgap=unit(7, "mm"),rightcols=c("HR",'ci','Pval', 'N_Incd_cases','N_Controls','N_Incd_cases_withVar','N_Controls_withVar'),
      rightlabs=c("HR",'95% CI','P', "Cases (N)", "Controls (N)", "Cases with CHIP (N)", "Controls with CHIP (N)"),plotwidth=unit(6.5, "cm"),
      col.square=c("grey", "grey", "red", "grey", "grey", "red"),
      #col.square=c("grey", "grey", "red", "grey", "grey", "red","grey", "grey", "red","grey", "grey", "red","grey", "grey", "red","grey", "grey", "red","grey", "grey", "red","grey", "grey", "red","grey", "grey", "red","grey", "grey", "red"),
      comb.random=F,print.Q=F,overall=F,hetstat=F,overall.hetstat=F,comb.fixed=F,print.byvar=F,addspace=T,digits.pval=3,xlab = "HR", digits=2)
dev.off()


###
UKB_Forest = read.table("/medpop/esp2/mzekavat/UKBB/PAD_CHIP/UKB_PAD_CHIP_assoc.txt",
                   header=T, as.is=T, stringsAsFactors=F, comment.char = '', sep="\t")
UKB_Forest$Cohort = "UK Biobank"

#Forest$Pval = as.character(Forest$Pr...z..)
Forest = UKB_Forest
Forest$Pval = Forest$Pr...z..


Forest$HR = Forest$exp.coef.
Forest$N_Controls = formatC(Forest$N_Controls, format = "f", big.mark = ",", drop0trailing = TRUE)
Forest$N_Incd_cases = formatC(Forest$N_Incd_cases, format = "f", big.mark = ",", drop0trailing = TRUE)
Forest$N_Incd_cases_withVar = formatC(Forest$N_Incd_cases_withVar, format = "f", big.mark = ",", drop0trailing = TRUE)
Forest$N_Controls_withVar = formatC(Forest$N_Controls_withVar, format = "f", big.mark = ",", drop0trailing = TRUE)
Forest <- Forest[order(Forest$Pr...z.. ),]


Forest$x = ifelse(Forest$x %in% c("CHIP", "NewCHIP"), "CHIP","Large CHIP")

library(meta)
library(grid)
library(scales)

Forest$Pval = formatC(Forest$Pval, format = "g", digits = 2)
Forest$x = ordered(Forest$x,levels=c("CHIP", "Large CHIP"))

Forest$Adjustment = ordered(Forest$Adjustment, levels=c("Unadjusted", "Sparsely Adjusted", "Fully Adjusted"))

Forest = Forest[order(Forest$x, Forest$Adjustment),]

summarymeta=metagen(coef,se.coef.,studlab=Adjustment,byvar = x,data=Forest,sm="HR")

pdf(paste("/medpop/esp2/mzekavat/UKBB/PAD_CHIP/Forest_CHIP.PAD.difftAdjustments.pdf",sep=""), width = 13, height= 3.5)
forest(summarymeta,leftcols =c("studlab"),leftlabs = c("          "),just='left',colgap.left=unit(0.1,"cm"),xlim=c(0.8, 6),at=c(0.8,1,2,4,6),
      smlab='Association of CHIP\nwith Incident Peripheral Artery Disease',colgap=unit(7, "mm"),rightcols=c("HR",'ci','Pval', 'N_Incd_cases','N_Controls','N_Incd_cases_withVar','N_Controls_withVar'),
      rightlabs=c("HR",'95% CI','P', "Cases (N)", "Controls (N)", "Cases with CHIP (N)", "Controls with CHIP (N)"),plotwidth=unit(6.5, "cm"),
      col.square=c("grey", "grey", "red", "grey", "grey", "red"),
      #col.square=c("grey", "grey", "red", "grey", "grey", "red","grey", "grey", "red","grey", "grey", "red","grey", "grey", "red","grey", "grey", "red","grey", "grey", "red","grey", "grey", "red","grey", "grey", "red","grey", "grey", "red"),
      comb.random=F,print.Q=F,overall=F,hetstat=F,overall.hetstat=F,comb.fixed=F,print.byvar=F,addspace=T,digits.pval=3,xlab = "HR", digits=2)
dev.off()

###########################################################################################################################################################################################################################################################################################################################
#--------- Propensity Score Analysis in UKB:

library(data.table)
library(survival)
library(tableone)

#########################################################################################################
### UKBB Phenos & QC:

pheno = fread("/medpop/esp2/mzekavat/UKBB/ukbb_PhenoFile.ALL_500k.UpdatedIncdPhenos_202020.txt")
pheno = data.frame(pheno)
pheno$IDs_toRemove_SampleQCv2 = ifelse((!(pheno$Submitted_Gender == pheno$Inferred_Gender) |  pheno$Second_deg_relOrHigher_toRemove ==1 | pheno$Non_Consented== 1),1,0) ## N:
pheno = pheno[which(is.na(pheno$IDs_toRemove_SampleQCv2 ) & pheno$inWES == 1),]
pheno = pheno[-which(pheno$Prev_Maryam_HemeCa_Phenos ==1),]
pheno = pheno[-which(pheno$Prev_Peripheral_vascular_disease ==1),]
pheno = pheno[which(!is.na(pheno$SmokingStatusv2) & !is.na(pheno$PC1) & !is.na(pheno$Townsend)),]

model_CHIP = glm(NewCHIP ~ age +age2 + Sex_numeric+SmokingStatusv2 + scale(Townsend)+ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, family=binomial(), data=pheno)
pheno$CHIP_propensityScore = predict(model_CHIP, type='response') 


library(ggplot2) #plotting system
library(grid)
library(scales)
library(ggrepel)
pheno$NewCHIP_label = ifelse(pheno$NewCHIP == 1, "+CHIP", ifelse(pheno$NewCHIP == 0, "-CHIP", NA) )
pheno$NewCHIP_label  = ordered(pheno$NewCHIP_label , levels=c("-CHIP", "+CHIP"))

p <- ggplot(pheno, aes(x=CHIP_propensityScore, fill=factor(NewCHIP_label)))
p <- p +geom_density(alpha=0.4, position='identity') + xlab('UKB CHIP Propensity Score')  +scale_fill_manual(breaks = c("-CHIP", "+CHIP"),values=c( "blue", "red"))
p <- p + theme(
	strip.background = element_blank(),
	strip.text.x = element_text(size = 14),
	strip.text.y = element_text(size = 14),
	axis.text.x = element_text(size=14, hjust=1, color = "black"),
	axis.text.y = element_text(size=14,hjust=1, color = "black"),
	axis.ticks =  element_line(colour = "black"), 
	axis.title.y= element_text(size=14),
	axis.title.x= element_text(size=14),
	legend.position="right", 
	legend.title = element_blank(),
	legend.text = element_text(size=14),
	panel.background = element_blank(), 
	panel.border = element_blank(), 
	panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank(), 
	panel.spacing = unit(1.0, "lines"), 
	plot.background = element_blank(), 
	plot.margin = unit(c(0.5,  0.5, 0.5, 0.5), "lines"),
	axis.line.x = element_line(colour = "black"),
	axis.line.y = element_line(colour = "black")
)+theme(legend.position="bottom")

pdf(paste("/medpop/esp2/mzekavat/Yale/CHIPpropScore.pdf",sep=""), width = 4, height= 4)
print(p)
dev.off()

model_LargeCHIP = glm(New_Large_CHIP_01 ~ age +age2 + Sex_numeric+SmokingStatusv2 + scale(Townsend)+ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, family=binomial(), data=pheno,na.action = na.exclude)
pheno$LargeCHIP_propensityScore = predict(model_LargeCHIP, type='response') 


pheno$NewLargeCHIP_label = ifelse(pheno$New_Large_CHIP_01 == 1, "+Large_CHIP", ifelse(pheno$New_Large_CHIP_01 == 0, "-Large_CHIP", NA) )
pheno$NewLargeCHIP_label  = ordered(pheno$NewLargeCHIP_label , levels=c("-Large_CHIP", "+Large_CHIP"))


p <- ggplot(pheno, aes(x=LargeCHIP_propensityScore, fill=factor(NewLargeCHIP_label)))
p <- p +geom_density(alpha=0.4, position='identity') + xlab('UKB Large CHIP Propensity Score')  +scale_fill_manual(breaks = c("-Large_CHIP", "+Large_CHIP"),values=c( "blue", "red"))
p <- p + theme(
	strip.background = element_blank(),
	strip.text.x = element_text(size = 14),
	strip.text.y = element_text(size = 14),
	axis.text.x = element_text(size=14, hjust=1, color = "black"),
	axis.text.y = element_text(size=14,hjust=1, color = "black"),
	axis.ticks =  element_line(colour = "black"), 
	axis.title.y= element_text(size=14),
	axis.title.x= element_text(size=14),
	legend.position="right", 
	legend.title = element_blank(),
	legend.text = element_text(size=14),
	panel.background = element_blank(), 
	panel.border = element_blank(), 
	panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank(), 
	panel.spacing = unit(1.0, "lines"), 
	plot.background = element_blank(), 
	plot.margin = unit(c(0.5,  0.5, 0.5, 0.5), "lines"),
	axis.line.x = element_line(colour = "black"),
	axis.line.y = element_line(colour = "black")
)+theme(legend.position="bottom")

pdf(paste("/medpop/esp2/mzekavat/Yale/LargeCHIPpropScore.pdf",sep=""), width = 4, height= 4)
print(p)
dev.off()

pheno$New_Large_CHIP = factor(pheno$New_Large_CHIP , levels=c("NO_CHIP",  "LARGE_CHIP","SMALL_CHIP"), ordered=FALSE)

CHIP_Names = c('NewCHIP','New_Large_CHIP')
pheno_list2 = c('Peripheral_vascular_disease')
vars = c(CHIP_Names)
summaryDF = data.frame()
for (i in 1:length(pheno_list2)){
	print(pheno_list2[i])
for (j in 1:length(vars)){
	tmp = pheno
	print(vars[j])
	#if(vars[j] %in%CHIP_Names){tmp = pheno}
	#if(vars[j] %in%c(CHUD_Names, CHUD_Normalized_Names)){tmp = pheno[-which(pheno$hasCHIP ==1),]}
indx = which(colnames(tmp) == vars[j])[1]

#removing prevalent cases
preval_col = which(colnames(tmp) == paste("Prev_",pheno_list2[i],sep=""))
incid_col = which(colnames(tmp) == paste("Incd_",pheno_list2[i],sep=""))

if (length(which(tmp[,preval_col] == 1))>0 ){tmp = tmp[-which(tmp[,preval_col] == 1),]}

if (length(which(tmp[,incid_col] ==1  &(tmp[,indx]>0 | tmp[,indx] == "LARGE_CHIP") & !is.na(tmp$age) & !is.na(tmp$PC1)  & !is.na(tmp$smok_detailed_))) > 4){

#running survival anal
tmp$SurvObj <- Surv(tmp[,which(colnames(tmp) == paste("FollowUp_",pheno_list2[i], sep=""))],  tmp[,which(colnames(tmp) == paste("Incd_",pheno_list2[i], sep=""))]== 1)

if(vars[j] == "NewCHIP"){df = as.data.frame(t(as.data.frame(summary( coxph(tmp$SurvObj ~ tmp[,indx]+scale(CHIP_propensityScore), data = tmp))$coeff[1,])))
}
else{df = as.data.frame(t(as.data.frame(summary( coxph(tmp$SurvObj ~ tmp[,indx]+scale(LargeCHIP_propensityScore), data = tmp))$coeff[1,])))
}
df$Adjustment = "PropensityScore"
df$y = pheno_list2[i]
df$x = vars[j]
df$N_Incd_cases = length(which(tmp[,incid_col] == 1 ))
df$N_Controls = length(which(tmp[,incid_col] == 0 ))
df$N_Incd_cases_withVar = length(which(tmp[,incid_col] == 1 & (tmp[,indx]>0 | tmp[,indx] == "LARGE_CHIP") ))
df$N_Controls_withVar = length(which(tmp[,incid_col] == 0 & (tmp[,indx]>0 | tmp[,indx] == "LARGE_CHIP") ))
summaryDF = rbind(summaryDF, df)

}
}
}


write.table(summaryDF,"/medpop/esp2/mzekavat/UKBB/PAD_CHIP/UKB_PAD_CHIP_assoc.PropScore_Adj.txt",
col.names = T, row.names = F, quote = F, sep = "\t")


CHIP_Prop = sum(pheno$NewCHIP) / length(which(!is.na(pheno$NewCHIP)))
LargeCHIP_Prop = sum(pheno$New_Large_CHIP_01,na.rm=TRUE) / length(which(!is.na(pheno$New_Large_CHIP_01)))

#-------CHIP weights for IPW:
weights_CHIP = rep(0, nrow(pheno))
id_chip1 = which(pheno$NewCHIP == 1)

weights_CHIP[id_chip1] = CHIP_Prop/pheno$CHIP_propensityScore[id_chip1]
weights_CHIP[-id_chip1] = (1-CHIP_Prop)/(1-pheno$CHIP_propensityScore[-id_chip1])
#--------Large CHIP weights for IPW:
weights_largeCHIP = rep(0, nrow(pheno))
id_largechip1 = which(pheno$New_Large_CHIP_01 == 1)

weights_largeCHIP[id_largechip1] = LargeCHIP_Prop/pheno$LargeCHIP_propensityScore[id_largechip1]
weights_largeCHIP[-id_largechip1] = (1-LargeCHIP_Prop)/(1-pheno$LargeCHIP_propensityScore[-id_largechip1])


CHIP_Names = c('NewCHIP','New_Large_CHIP')
pheno_list2 = c('Peripheral_vascular_disease')
vars = c(CHIP_Names)
summaryDF = data.frame()
for (i in 1:length(pheno_list2)){
	print(pheno_list2[i])
for (j in 1:length(vars)){
	tmp = pheno
	print(vars[j])
	#if(vars[j] %in%CHIP_Names){tmp = pheno}
	#if(vars[j] %in%c(CHUD_Names, CHUD_Normalized_Names)){tmp = pheno[-which(pheno$hasCHIP ==1),]}
indx = which(colnames(tmp) == vars[j])[1]

#removing prevalent cases
preval_col = which(colnames(tmp) == paste("Prev_",pheno_list2[i],sep=""))
incid_col = which(colnames(tmp) == paste("Incd_",pheno_list2[i],sep=""))

if (length(which(tmp[,preval_col] == 1))>0 ){tmp = tmp[-which(tmp[,preval_col] == 1),]}

if (length(which(tmp[,incid_col] ==1  &(tmp[,indx]>0 | tmp[,indx] == "LARGE_CHIP") & !is.na(tmp$age) & !is.na(tmp$PC1)  & !is.na(tmp$smok_detailed_))) > 4){

#running survival anal
tmp$SurvObj <- Surv(tmp[,which(colnames(tmp) == paste("FollowUp_",pheno_list2[i], sep=""))],  tmp[,which(colnames(tmp) == paste("Incd_",pheno_list2[i], sep=""))]== 1)

if(vars[j] == "NewCHIP"){df = as.data.frame(t(as.data.frame(summary( coxph(tmp$SurvObj ~ tmp[,indx], data = tmp, weights=weights_CHIP))$coeff[1,])))
}
else{df = as.data.frame(t(as.data.frame(summary( coxph(tmp$SurvObj ~ tmp[,indx], data = tmp, weights=weights_largeCHIP))$coeff[1,])))
}
df$Adjustment = "IPTW"
df$y = pheno_list2[i]
df$x = vars[j]
df$N_Incd_cases = length(which(tmp[,incid_col] == 1 ))
df$N_Controls = length(which(tmp[,incid_col] == 0 ))
df$N_Incd_cases_withVar = length(which(tmp[,incid_col] == 1 & (tmp[,indx]>0 | tmp[,indx] == "LARGE_CHIP") ))
df$N_Controls_withVar = length(which(tmp[,incid_col] == 0 & (tmp[,indx]>0 | tmp[,indx] == "LARGE_CHIP") ))
summaryDF = rbind(summaryDF, df)

}
}
}

write.table(summaryDF,"/medpop/esp2/mzekavat/UKBB/PAD_CHIP/UKB_PAD_CHIP_assoc.IPTW.txt",
col.names = T, row.names = F, quote = F, sep = "\t")


UKB_Forest1 = read.table("/medpop/esp2/mzekavat/UKBB/PAD_CHIP/UKB_PAD_CHIP_assoc.PropScore_Adj.txt",
                   header=T, as.is=T, stringsAsFactors=F, comment.char = '', sep="\t")
UKB_Forest1$Cohort = "UK Biobank"

UKB_Forest2 = read.table("/medpop/esp2/mzekavat/UKBB/PAD_CHIP/UKB_PAD_CHIP_assoc.IPTW.txt",
                   header=T, as.is=T, stringsAsFactors=F, comment.char = '', sep="\t")
UKB_Forest2$Cohort = "UK Biobank"

#Forest$Pval = as.character(Forest$Pr...z..)
Forest = rbind(UKB_Forest1, UKB_Forest2)
Forest$Pval = Forest$Pr...z..


Forest$HR = Forest$exp.coef.
Forest$N_Controls = formatC(Forest$N_Controls, format = "f", big.mark = ",", drop0trailing = TRUE)
Forest$N_Incd_cases = formatC(Forest$N_Incd_cases, format = "f", big.mark = ",", drop0trailing = TRUE)
Forest$N_Incd_cases_withVar = formatC(Forest$N_Incd_cases_withVar, format = "f", big.mark = ",", drop0trailing = TRUE)
Forest$N_Controls_withVar = formatC(Forest$N_Controls_withVar, format = "f", big.mark = ",", drop0trailing = TRUE)
Forest <- Forest[order(Forest$Pr...z.. ),]


Forest$x = ifelse(Forest$x %in% c("CHIP", "NewCHIP"), "CHIP","Large CHIP")

library(meta)
library(grid)
library(scales)

Forest$Pval = formatC(Forest$Pval, format = "g", digits = 2)
Forest$x = ordered(Forest$x,levels=c("CHIP", "Large CHIP"))
Forest$Adjustment = ifelse(Forest$Adjustment == "IPTW", "IPTW", "Propensity Score Adjustment")
Forest$Adjustment = ordered(Forest$Adjustment, levels=c("Propensity Score Adjustment", "IPTW"))

Forest = Forest[order(Forest$x, Forest$Adjustment),]

summarymeta=metagen(coef,se.coef.,studlab=Adjustment,byvar = x,data=Forest,sm="HR")

pdf(paste("/medpop/esp2/mzekavat/UKBB/PAD_CHIP/Forest_CHIP.propScore_Adjustments.pdf",sep=""), width = 11, height= 3)
forest(summarymeta,leftcols =c("studlab"),leftlabs = c("          "),just='left',colgap.left=unit(0.1,"cm"),xlim=c(0.8, 6),at=c(0.8,1,2,4,6),
      smlab='Association of CHIP\nwith Incident Peripheral Artery Disease',colgap=unit(7, "mm"),rightcols=c("HR",'ci','Pval'),
      rightlabs=c("HR",'95% CI','P'),plotwidth=unit(6.5, "cm"),
      col.square=c("grey", "grey", "grey", "grey"),
      #col.square=c("grey", "grey", "red", "grey", "grey", "red","grey", "grey", "red","grey", "grey", "red","grey", "grey", "red","grey", "grey", "red","grey", "grey", "red","grey", "grey", "red","grey", "grey", "red","grey", "grey", "red"),
      comb.random=F,print.Q=F,overall=F,hetstat=F,overall.hetstat=F,comb.fixed=F,print.byvar=F,addspace=T,digits.pval=3,xlab = "HR", digits=2)
dev.off()


