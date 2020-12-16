#Mendelian randomization analysis (CHIP on PAD)
#Maryam Zekavat
#Dec 15, 2020

###---------------------------------------------reading in CHIP gwas and PAD gwas; identifying variants with CHIP P<0.001 in common between the CHIP and PAD GWAS:

gwas <- fread("/medpop/esp2/mesbah/Meta_GWAS/CHIP_75kTOPMED_48kUKB.tsv.gz", header=TRUE, sep="\t")
gwas = data.frame(gwas)
gwas$P = 10^(-gwas$neg_log_pvalue)

PAD = fread("/medpop/esp2/dklarin/MVP/PAD/CLEANED/anno/CLEANED.MVP.EUR.PAD.results.anno.csv", header=TRUE)
PAD = data.frame(PAD)
gwas = gwas[which(gwas$P<=0.001),]

gwas = gwas[which(gwas$rsid %in% PAD$ID),]

write.table(gwas,"/medpop/esp2/mzekavat/Miscellaneous_Tables/GWAS_Results/CHIP_GWAS/logreg_wald_hasCHIP.p_.001.inMVP.tsv",
col.names = T, row.names = F, quote = F, sep = "\t")

###---------------------------------------------Using Plink-1.9 to perform ld-clumping prioritized by p-value using the 1000G-european reference panel:
use .plink-1.90b

plink --bfile /medpop/esp2/akhil/topmed/lipid_mashr/1kgEUR/1000G.EUR \
    --maf 0.00005 \
    --clump /medpop/esp2/mzekavat/Miscellaneous_Tables/GWAS_Results/CHIP_GWAS/logreg_wald_hasCHIP.p_.001.inMVP.tsv \
    --clump-p1 0.001 \
    --clump-p2 1 \
    --clump-r2 0.1 \
    --clump-kb 1000 \
    --clump-snp-field rsid \
    --out /medpop/esp2/mzekavat/Miscellaneous_Tables/GWAS_Results/CHIP_GWAS/logreg_wald_hasCHIP.p_.001.inMVP.Clumped #--clump: 427 clumps formed from 12661 top variants

###--------------------------------------------- Now merging list of significant, independent SNPs with CHIP GWAS summary statistics and PAD GWAS summary statistics:
CHIP <- fread("/medpop/esp2/mzekavat/Miscellaneous_Tables/GWAS_Results/CHIP_GWAS/logreg_wald_hasCHIP.p_.001.inMVP.Clumped.clumped", header=TRUE)
CHIP = data.frame(CHIP)

CHIPgwas = fread("/medpop/esp2/mzekavat/Miscellaneous_Tables/GWAS_Results/CHIP_GWAS/logreg_wald_hasCHIP.p_.001.inMVP.tsv", header=TRUE, sep="\t")
CHIPgwas  = data.frame(CHIPgwas )

CHIPgwas = CHIPgwas[which(CHIPgwas$rsid %in% CHIP$SNP),]

PAD = fread("/medpop/esp2/dklarin/MVP/PAD/CLEANED/anno/CLEANED.MVP.EUR.PAD.results.anno.csv", header=TRUE)
PAD = data.frame(PAD)


colnames(CHIPgwas) = paste("CHIP_", colnames(CHIPgwas), sep="")
colnames(PAD) = paste("PAD_",colnames(PAD),sep="")

merged = merge(CHIPgwas,PAD,by.x=c(3), by.y=c(16), all.x=TRUE)

merged$new_PAD_beta = ifelse(merged$CHIP_alt ==merged$PAD_EFFECT_ALLELE, merged$PAD_BETA,
                    ifelse(merged$CHIP_alt ==merged$PAD_OTHER_ALLELE, -1*merged$PAD_BETA, NA))


write.table(merged,"/medpop/esp2/mzekavat/Miscellaneous_Tables/GWAS_Results/CHIP_GWAS/PAD_CHIP.MRanal.txt",
col.names = T, row.names = F, quote = F, sep = "\t")

###---------------------------------------------removing variants associated with Smoking Initiation and filtering to CHIP-P<0.0001:
merged = fread("/medpop/esp2/mzekavat/Miscellaneous_Tables/GWAS_Results/CHIP_GWAS/PAD_CHIP.MRanal.txt", header=TRUE, sep="\t")
merged  = data.frame(merged )

Smokinggwas <- fread("/medpop/esp2/mzekavat/Miscellaneous_Tables/GWAS_Results/Smoking/SmokingInitiation.WithoutUKB.txt.gz", header=TRUE, sep="\t")
Smokinggwas = data.frame(Smokinggwas)
colnames(Smokinggwas)= paste('Smoking_', colnames(Smokinggwas), sep="")

mergedv2 = merge(merged, Smokinggwas, by.x=c(1), by.y=c(3), all.x=TRUE)

write.table(mergedv2,"/medpop/esp2/mzekavat/Miscellaneous_Tables/GWAS_Results/CHIP_GWAS/PAD_CHIP.MRanal.withSmokingAnnot.txt",
col.names = T, row.names = F, quote = F, sep = "\t") #40var / 3160 var

mergedv2 = mergedv2[-which(mergedv2$Smoking_PVALUE<0.05),]
mergedv2 = mergedv2[which(mergedv2$CHIP_P<1e-4),]

write.table(mergedv2,"/medpop/esp2/mzekavat/Miscellaneous_Tables/GWAS_Results/CHIP_GWAS/PAD_CHIP.MRanal.withSmokingAnnot.minusSmokingSignVar.CHIPp0.0001.txt",
col.names = T, row.names = F, quote = F, sep = "\t") #40var / 3160 var

###---------------------------------------------Performing MR:

merged = fread("/medpop/esp2/mzekavat/Miscellaneous_Tables/GWAS_Results/CHIP_GWAS/PAD_CHIP.MRanal.withSmokingAnnot.minusSmokingSignVar.CHIPp0.0001.txt", header=TRUE, sep="\t")
merged  = data.frame(merged )

library(MendelianRandomization)
library(data.table)
mergedv2 = fread("/medpop/esp2/mzekavat/Miscellaneous_Tables/GWAS_Results/CHIP_GWAS/PAD_CHIP.MRanal.withSmokingAnnot.minusSmokingSignVar.CHIPp0.0001.txt")
mergedv2 = data.frame(mergedv2)
merged = mergedv2[which(!is.na(mergedv2$new_PAD_beta)),]
mergedv2 = merged[which(merged$CHIP_P<1e-4),] #404 variants
mergedv2 = mergedv2[!duplicated(mergedv2$CHIP_rsid),]

###---------------------------------------------Making MR scatter plot (Fig 6c)
p=mr_plot(
  mr_input(bx = mergedv2$CHIP_beta, bxse = mergedv2$CHIP_stderr_beta,by = mergedv2$new_PAD_beta,byse =  mergedv2$PAD_SE),
  error = TRUE,
  line = "ivw",
  orientate = FALSE,
  interactive = FALSE,
  labels = FALSE
)
pdf(paste("/medpop/esp2/mzekavat/Miscellaneous_Tables/GWAS_Results/CHIP_GWAS/CHIP_PAD_scatplot.pdf",sep=""), width = 6, height= 5)
print(p)
dev.off()

snpList=mergedv2$CHIP_rsid[which(mergedv2$CHIP_beta*mergedv2$new_PAD_beta>0 & mergedv2$PAD_PVAL<0.05)]
res=phenoscanner(snpquery=snpList,pvalue = 1e-05, proxies = "None", r2 = 0.8)

write.table(res$results,"/medpop/esp2/mzekavat/Miscellaneous_Tables/GWAS_Results/CHIP_GWAS/PAD_CHIP.phenoScannerVar_withSimilarEffects.txt",
col.names = T, row.names = F, quote = F, sep = "\t") #40var / 3160 var

###---------------------------------------------Performing the contamination mixture MR model:

mr_conmix(mr_input(bx = mergedv2$CHIP_beta, bxse = mergedv2$CHIP_stderr_beta,by = mergedv2$new_PAD_beta,byse =  mergedv2$PAD_SE),psi = 3, CIMin = -1, CIMax = 1, CIStep = 0.001)

Contamination mixture method
(Standard deviation of invalid estimands = 3)
Number of Variants : 404 
------------------------------------------------------------------
 Method Estimate 95% CI       p-value
 ConMix    0.017 0.003, 0.032  0.0244
-------------------------------------------------------

###---------------------------------------------Performing other MR models:
mr_allmethods(mr_input(bx = mergedv2$CHIP_beta, bxse = mergedv2$CHIP_stderr_beta,by = mergedv2$new_PAD_beta,byse =  mergedv2$PAD_SE),method = "all")
> mr_allmethods(mr_input(bx = mergedv2$CHIP_beta, bxse = mergedv2$CHIP_stderr_beta,by = mergedv2$new_PAD_beta,byse =  mergedv2$PAD_SE),method = "all")
                    Method Estimate Std Error 95% CI        P-value
             Simple median    0.023     0.012   0.000 0.046   0.049
           Weighted median    0.019     0.012  -0.005 0.043   0.121
 Penalized weighted median    0.015     0.012  -0.009 0.038   0.229
                                                                   
                       IVW    0.018     0.008   0.002 0.033   0.028
             Penalized IVW    0.017     0.008   0.001 0.032   0.034
                Robust IVW    0.016     0.008   0.000 0.032   0.053
      Penalized robust IVW    0.016     0.008   0.000 0.032   0.051
                                                                   
                  MR-Egger   -0.010     0.014  -0.038 0.018   0.480
               (intercept)    0.004     0.002   0.001 0.008   0.018
        Penalized MR-Egger   -0.011     0.014  -0.038 0.016   0.406
               (intercept)    0.005     0.002   0.001 0.008   0.013
           Robust MR-Egger   -0.011     0.013  -0.036 0.014   0.387
               (intercept)    0.005     0.002   0.001 0.008   0.019
 Penalized robust MR-Egger   -0.012     0.013  -0.037 0.013   0.363
               (intercept)    0.005     0.002   0.001 0.008   0.016
> 

###--------------------------------------------- Performing MR-lasso model:

mr_lasso(mr_input(bx = mergedv2$CHIP_beta, bxse = mergedv2$CHIP_stderr_beta,by = mergedv2$new_PAD_beta,byse =  mergedv2$PAD_SE))
MR-Lasso method 

Number of variants : 404 
Number of valid instruments : 398 
Tuning parameter : 0.1323845 
------------------------------------------------------------------
 Exposure Estimate Std Error 95% CI       p-value
 exposure    0.016     0.008 0.001, 0.031   0.043
------------------------------------------------------------------

###--------------------------------------------- Performing MR-RAPS:

library(mr.raps)
res=mr.raps.mle.all(mergedv2$CHIP_beta, mergedv2$new_PAD_beta, mergedv2$CHIP_stderr_beta, mergedv2$PAD_SE)
res$beta.z=res$beta.hat / res$beta.se
res$beta.p=2*pnorm(-abs(res$beta.z))

res
  over.dispersion loss.function   beta.hat     beta.se   beta.z      beta.p
1           FALSE            l2 0.02291196 0.007932661 2.888306 0.003873225
2           FALSE         huber 0.02051250 0.008136468 2.521057 0.011700277
3           FALSE         tukey 0.02062627 0.008136585 2.535004 0.011244607
4            TRUE            l2 0.02186268 0.008160757 2.679001 0.007384215
5            TRUE         huber 0.02051253 0.008137013 2.520892 0.011705764
6            TRUE         tukey 0.02062637 0.008137156 2.534838 0.011249930

###--------------------------------------------- Making forest plot of MR associations:
str(p$data)
'data.frame':   4 obs. of  4 variables:
 $ snps     : Factor w/ 4 levels "Contamination mixture estimate",..: 4 3 2 1
 $ estimates: num  0.0175 0.0188 0.0179 0.022
 $ CI_lower : num  0.00186 -0.00497 0.0018 0.00202
 $ CI_upper : num  0.0332 0.0426 0.0341 0.0

p=mr_forest(mr_input(bx = mergedv2$CHIP_beta, bxse = mergedv2$CHIP_stderr_beta,by = mergedv2$new_PAD_beta,byse =  mergedv2$PAD_SE), methods = c("ivw", "wmedian", 'lasso', 'maxlik','conmix'), snp_estimates = FALSE)

df = data.frame(snps = c("MR-lasso", "MR-RAPS (l2-loss function)"), estimates = c(0.016, 0.0229), CI_lower = c(0.001,0.00735),CI_upper =c(0.031,0.0384))
p$data = rbind(p$data, df)
p$data$snps = ordered(p$data$snps, levels=c('Contamination mixture estimate','MR-RAPS (l2-loss function)','Maximum likelihood estimate','MR-lasso','Weighted median estimate','IVW estimate'))
pdf(paste("/medpop/esp2/mzekavat/Miscellaneous_Tables/GWAS_Results/CHIP_GWAS/CHIP_PAD_forestplot.v2.pdf",sep=""), width = 6, height= 3)
print(p)
dev.off()

###--------------------------------------------- Sensitivity analyses using TwoSampleMR package:

library(TwoSampleMR)
outcome_dat <- read_outcome_data(
    filename = "/medpop/esp2/mzekavat/Miscellaneous_Tables/GWAS_Results/CHIP_GWAS/PAD_CHIP.MRanal.withSmokingAnnot.minusSmokingSignVar.CHIPp0.0001.txt",
    sep = "\t",
    snp_col = "CHIP_rsid",
    beta_col = "new_PAD_beta",
    se_col = "PAD_SE",
    effect_allele_col = "CHIP_alt",
    other_allele_col = "CHIP_ref",
    eaf_col = "CHIP_alt_allele_freq",
    pval_col = "PAD_PVAL"

    )
#    eaf_col = "a1_freq",

exposure_dat <- read_exposure_data(
    filename = "/medpop/esp2/mzekavat/Miscellaneous_Tables/GWAS_Results/CHIP_GWAS/PAD_CHIP.MRanal.withSmokingAnnot.minusSmokingSignVar.CHIPp0.0001.txt",
    sep = "\t",
    snp_col = "CHIP_rsid",
    beta_col = "CHIP_beta",
    se_col = "CHIP_stderr_beta",
    effect_allele_col = "CHIP_alt",
    other_allele_col = "CHIP_ref",
    eaf_col = "CHIP_alt_allele_freq",
    pval_col = "CHIP_P"
    )

dat = harmonise_data(outcome_dat, exposure_dat)

### Performing MR using mr_raps and mr_ivw and making scatter plots:
res <- mr(dat, method_list=c("mr_ivw", "mr_simple_median", "mr_weighted_median"))


library(ggplot2) #plotting system
library(grid)
library(reshape2)
library(scales)
library(ggrepel)
library(plotly)
library(meta)
library(tidyr)
library(ggsci)
library(data.table)
p1 <- mr_scatter_plot(res, dat)#+theme(axis.text.y=element_text(size=18, hjust=1, color='black'),axis.text.x=element_text(size=18, hjust=1, color='black'),axis.title.x= element_text(size=18))+ xlab("SNP effect on CHIP")+ ylab("SNP effect on PAD")

pdf(paste("/medpop/esp2/mzekavat/Miscellaneous_Tables/GWAS_Results/CHIP_GWAS/CHIP_PAD_scatterplot.pdf",sep=""), width = 4, height= 5)
p1[[1]]
dev.off()

###--------------------------------------------- Performing Sensitivity Analyses:



#--------------------------------------------- Evaluating MR Assumptions using UKB: 

"CHIP_p1e4_PRS"











