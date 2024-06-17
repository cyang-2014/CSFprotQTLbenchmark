rm(list = ls())
library(dplyr)
library(data.table)
library(reshape2)
library(openxlsx)
library(Biobase)
library(ggplot2)
library(cowplot)

library(foreach)

#######################################
## only test the proteins with sentinel variants from a union set I curated on 2/28/2024
DIRout_prep <- '/03-DryLab/04-Analyses/2021_Multi-Omics_CC/2022_MultiTissue-MultiOmics_Chengran/t16_benchmarkQTL/run00_genoCSFlongs/'
list.files(DIRout_prep)

unionSets.dt.uniqID.protFound <- fread(paste0(DIRout_prep, 'f06a_union_keepUniqueProtVarBlock_pairs.csv'))
unionSets.dt.uniqID.protFound
uniqueN(unionSets.dt.uniqID.protFound$protVarBlockID)
uniqueN(unionSets.dt.uniqID.protFound$Analytes)

table(unionSets.dt.uniqID.protFound$setName)
table(unionSets.dt.uniqID.protFound[chr.var == 'chr3']$LDname, unionSets.dt.uniqID.protFound[chr.var == 'chr3']$setName) # LDname == '394_START' from DW2022 as expected

#######################################
## set up the plink2 glm command
genodir <- '/03-DryLab/04-Analyses/2021_Multi-Omics_CC/2022_MultiTissue-MultiOmics_Chengran/t16_benchmarkQTL/run00_genoCSFlongs/st01_extractFEB2024/'
list.files(genodir)
filename <- paste0(genodir, 'CSF_pQTL_benchmark_FEB2024_arrayOptimized_final')

phenodir <- '/03-DryLab/04-Analyses/2021_Multi-Omics_CC/2022_MultiTissue-MultiOmics_Chengran/t16_benchmarkQTL/run01_CSFlongsLOG10z/d00_input/'
list.files(phenodir)
pheno <- paste0(phenodir, 'f02a_EURgeno_pheno_2249pariticipants_7008somamers.txt')

covardir <- '/03-DryLab/04-Analyses/2021_Multi-Omics_CC/2022_MultiTissue-MultiOmics_Chengran/t16_benchmarkQTL/run01_CSFlongsLOG10z/d00_input/'
covarFILE<- paste0(covardir, 'f04_EURgeno_covar_2249pariticipants_10gPC60peer.txt')


outdir <- '/03-DryLab/04-Analyses/2021_Multi-Omics_CC/2022_MultiTissue-MultiOmics_Chengran/t16_benchmarkQTL/run01_CSFlongsLOG10z/d01_plink2OUTwiPEER/t02_sentinelVars/output/'
list.files(outdir)

## loop each protein-variant pair for plink2 glm test
# for(i in 1:10) {
for(i in 1:nrow(unionSets.dt.uniqID.protFound)) {
  ## Run PLINK for each analyte using the covariates from argument and any additional arguments
  cmd_plink2 <- paste0('/home/yangc/myprogram/plink2_20230607 --bfile ', filename,
                       ' --glm log10 hide-covar',
                       ' --covar ', covarFILE,
                       ' --covar-variance-standardize',
                       ' --pheno ', pheno,
                       ' --pheno-name ', unionSets.dt.uniqID.protFound$Analytes[i],
                       ' --snp ', unionSets.dt.uniqID.protFound$snpID[i],
                       ' --out ', outdir, 'T01_CSFplink2_sentinelVar', 
                       unionSets.dt.uniqID.protFound$protVarBlockID[i],
                       ' --threads 24')
  system(cmd_plink2)
  cat(i, 'finished \n')

}



### to assemble into one table for the downstream comparisons between models
# plink2sentOUT.all <- foreach(i = 1:10, .combine = 'rbind') %do% {
# # plink2sentOUT.all <- foreach(i = 1:nrow(unionSets.dt.uniqID.protFound), .combine = 'rbind') %do% {
#   plink2sentOUT <- fread(paste0(outdir, 'T01_CSFplink2_sentinelVar',
#                                 unionSets.dt.uniqID.protFound$protVarBlockID[i], '.',
#                                 unionSets.dt.uniqID.protFound$Analytes[i], '.glm.linear'))
#   plink2sentOUT$protID.new <- unionSets.dt.uniqID.protFound$Analytes[i]
#   return(plink2sentOUT)
# }
# 
# 
# plink2sentOUT.all





