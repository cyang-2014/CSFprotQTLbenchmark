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
## load in the plink2 glm output


### T01
outdir.T01 <- '/03-DryLab/04-Analyses/2021_Multi-Omics_CC/2022_MultiTissue-MultiOmics_Chengran/t16_benchmarkQTL/run01_CSFlongsLOG10z/d01_plink2OUTwiPEER/t02_sentinelVars/output/'
list.files(outdir.T01)

## to assemble into one table for the downstream comparisons between models
# plink2sentOUT.all_T01 <- foreach(i = 1:10, .combine = 'rbind') %do% {
plink2sentOUT.all_T01 <- foreach(i = 1:nrow(unionSets.dt.uniqID.protFound), .combine = 'rbind') %do% {
  plink2sentOUT <- fread(paste0(outdir.T01, 'T01_CSFplink2_sentinelVar',
                                unionSets.dt.uniqID.protFound$protVarBlockID[i], '.',
                                unionSets.dt.uniqID.protFound$Analytes[i], '.glm.linear'))
  plink2sentOUT$protID.new <- unionSets.dt.uniqID.protFound$Analytes[i]
  plink2sentOUT$protVarBlockID <- unionSets.dt.uniqID.protFound$protVarBlockID[i]
  return(plink2sentOUT)
}


plink2sentOUT.all_T01
plink2sentOUT.all_T01.compare <- merge(plink2sentOUT.all_T01, 
                                       unionSets.dt.uniqID.protFound[, list(protVarBlockID, setName, LDname, LP)], 
                                       by = 'protVarBlockID')
plink2sentOUT.all_T01.compare[LDname == '394_START']
cor.test(plink2sentOUT.all_T01.compare[LDname == '394_START']$LP,
         plink2sentOUT.all_T01.compare[LDname == '394_START']$LOG10_P)

ggplot(plink2sentOUT.all_T01.compare[LDname == '394_START'], aes(LP, LOG10_P)) + geom_point()
ggplot(plink2sentOUT.all_T01.compare[LDname == '394_START'], aes(LP, LOG10_P)) + 
  geom_point() +
  xlim(0, 20) +
  ylim(0, 20)



### T02
outdir.T02 <- '/03-DryLab/04-Analyses/2021_Multi-Omics_CC/2022_MultiTissue-MultiOmics_Chengran/t16_benchmarkQTL/run01_CSFlongsLOG10z/d02_plink2OUTnoPEER/t02_sentinelVars/output/'
list.files(outdir.T02)

## to assemble into one table for the downstream comparisons between models
# plink2sentOUT.all_T02 <- foreach(i = 1:10, .combine = 'rbind') %do% {
plink2sentOUT.all_T02 <- foreach(i = 1:nrow(unionSets.dt.uniqID.protFound), .combine = 'rbind') %do% {
  plink2sentOUT <- fread(paste0(outdir.T02, 'T02_CSFplink2_sentinelVar',
                                unionSets.dt.uniqID.protFound$protVarBlockID[i], '.',
                                unionSets.dt.uniqID.protFound$Analytes[i], '.glm.linear'))
  plink2sentOUT$protID.new <- unionSets.dt.uniqID.protFound$Analytes[i]
  plink2sentOUT$protVarBlockID <- unionSets.dt.uniqID.protFound$protVarBlockID[i]
  return(plink2sentOUT)
}


plink2sentOUT.all_T02


plink2sentOUT.all_T02.compare <- merge(plink2sentOUT.all_T02, 
                                       unionSets.dt.uniqID.protFound[, list(protVarBlockID, setName, LDname, LP)], 
                                       by = 'protVarBlockID')
plink2sentOUT.all_T02.compare[LDname == '394_START']
table(plink2sentOUT.all_T02.compare[LDname == '394_START']$setName)
cor.test(plink2sentOUT.all_T02.compare[LDname == '394_START']$LP,
         plink2sentOUT.all_T02.compare[LDname == '394_START']$LOG10_P)

ggplot(plink2sentOUT.all_T02.compare[LDname == '394_START'], aes(LP, LOG10_P)) + geom_point()
ggplot(plink2sentOUT.all_T02.compare[LDname == '394_START'], aes(LP, LOG10_P)) + 
  geom_point() +
  xlim(0, 20) +
  ylim(0, 20)

## maybe PPMI in the meta-DW2022 raise these chr3-LD394 trans QTL pass the study-wide significance?


### T05
outdir.T05 <- '/03-DryLab/04-Analyses/2021_Multi-Omics_CC/2022_MultiTissue-MultiOmics_Chengran/t16_benchmarkQTL/run02_CSFlongsINVR/d05_plink2OUTwiPEER/t02_sentinelVars/output/'
list.files(outdir.T05)

## to assemble into one table for the downstream comparisons between models
# plink2sentOUT.all_T05 <- foreach(i = 1:10, .combine = 'rbind') %do% {
plink2sentOUT.all_T05 <- foreach(i = 1:nrow(unionSets.dt.uniqID.protFound), .combine = 'rbind') %do% {
  plink2sentOUT <- fread(paste0(outdir.T05, 'T05_CSFplink2_sentinelVar',
                                unionSets.dt.uniqID.protFound$protVarBlockID[i], '.',
                                unionSets.dt.uniqID.protFound$Analytes[i], '.glm.linear'))
  plink2sentOUT$protID.new <- unionSets.dt.uniqID.protFound$Analytes[i]
  plink2sentOUT$protVarBlockID <- unionSets.dt.uniqID.protFound$protVarBlockID[i]
  return(plink2sentOUT)
}


plink2sentOUT.all_T05
plink2sentOUT.all_T05.compare <- merge(plink2sentOUT.all_T05, 
                                       unionSets.dt.uniqID.protFound[, list(protVarBlockID, setName, LDname, LP)], 
                                       by = 'protVarBlockID')
plink2sentOUT.all_T05.compare[LDname == '394_START']
cor.test(plink2sentOUT.all_T05.compare[LDname == '394_START']$LP,
         plink2sentOUT.all_T05.compare[LDname == '394_START']$LOG10_P)

ggplot(plink2sentOUT.all_T05.compare[LDname == '394_START'], aes(LP, LOG10_P)) + geom_point()
ggplot(plink2sentOUT.all_T05.compare[LDname == '394_START'], aes(LP, LOG10_P)) + 
  geom_point() +
  xlim(0, 20) +
  ylim(0, 20)


### T06
outdir.T06 <- '/03-DryLab/04-Analyses/2021_Multi-Omics_CC/2022_MultiTissue-MultiOmics_Chengran/t16_benchmarkQTL/run02_CSFlongsINVR/d06_plink2OUTnoPEER/t02_sentinelVars/output/'
list.files(outdir.T06)

## to assemble into one table for the downstream comparisons between models
# plink2sentOUT.all_T06 <- foreach(i = 1:10, .combine = 'rbind') %do% {
plink2sentOUT.all_T06 <- foreach(i = 1:nrow(unionSets.dt.uniqID.protFound), .combine = 'rbind') %do% {
  plink2sentOUT <- fread(paste0(outdir.T06, 'T06_CSFplink2_sentinelVar',
                                unionSets.dt.uniqID.protFound$protVarBlockID[i], '.',
                                unionSets.dt.uniqID.protFound$Analytes[i], '.glm.linear'))
  plink2sentOUT$protID.new <- unionSets.dt.uniqID.protFound$Analytes[i]
  plink2sentOUT$protVarBlockID <- unionSets.dt.uniqID.protFound$protVarBlockID[i]
  return(plink2sentOUT)
}


plink2sentOUT.all_T06
plink2sentOUT.all_T06.compare <- merge(plink2sentOUT.all_T06, 
                                       unionSets.dt.uniqID.protFound[, list(protVarBlockID, setName, LDname, LP)], 
                                       by = 'protVarBlockID')
plink2sentOUT.all_T06.compare[LDname == '394_START']
cor.test(plink2sentOUT.all_T06.compare[LDname == '394_START']$LP,
         plink2sentOUT.all_T06.compare[LDname == '394_START']$LOG10_P)

ggplot(plink2sentOUT.all_T06.compare[LDname == '394_START'], aes(LP, LOG10_P)) + geom_point()
ggplot(plink2sentOUT.all_T06.compare[LDname == '394_START'], aes(LP, LOG10_P)) + 
  geom_point() +
  xlim(0, 20) +
  ylim(0, 20)


### write out the assembled file from each model
new.outdir.T01 <- '/03-DryLab/04-Analyses/2021_Multi-Omics_CC/2022_MultiTissue-MultiOmics_Chengran/t16_benchmarkQTL/run01_CSFlongsLOG10z/d01_plink2OUTwiPEER/t02_sentinelVars/'
list.files(new.outdir.T01)

fwrite(plink2sentOUT.all_T01, paste0(new.outdir.T01, 'assembled_out_T01.csv'))


new.outdir.T02 <- '/03-DryLab/04-Analyses/2021_Multi-Omics_CC/2022_MultiTissue-MultiOmics_Chengran/t16_benchmarkQTL/run01_CSFlongsLOG10z/d02_plink2OUTnoPEER/t02_sentinelVars/'
list.files(new.outdir.T02)

fwrite(plink2sentOUT.all_T02, paste0(new.outdir.T02, 'assembled_out_T02.csv'))

new.outdir.T05 <- '/03-DryLab/04-Analyses/2021_Multi-Omics_CC/2022_MultiTissue-MultiOmics_Chengran/t16_benchmarkQTL/run02_CSFlongsINVR/d05_plink2OUTwiPEER/t02_sentinelVars/'
list.files(new.outdir.T05)

fwrite(plink2sentOUT.all_T05, paste0(new.outdir.T05, 'assembled_out_T05.csv'))

new.outdir.T06 <- '/03-DryLab/04-Analyses/2021_Multi-Omics_CC/2022_MultiTissue-MultiOmics_Chengran/t16_benchmarkQTL/run02_CSFlongsINVR/d06_plink2OUTnoPEER/t02_sentinelVars/'
list.files(new.outdir.T06)

fwrite(plink2sentOUT.all_T06, paste0(new.outdir.T06, 'assembled_out_T06.csv'))

























