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
## load in the regenie output

### T03
outdir.T03 <- '/03-DryLab/04-Analyses/2021_Multi-Omics_CC/2022_MultiTissue-MultiOmics_Chengran/t16_benchmarkQTL/run01_CSFlongsLOG10z/d03_regenieOUTwiPEER/t02_sentinelVars/outputS2/'
list.files(outdir.T03)
## to assemble into one table for the downstream comparisons between models
# regeniesentOUT.all_T03 <- foreach(i = 1:10, .combine = 'rbind') %do% {
regeniesentOUT.all_T03 <- foreach(i = 1:nrow(unionSets.dt.uniqID.protFound), .combine = 'rbind') %do% {
  regeniesentOUT <- fread(paste0(outdir.T03, 'T03_WU_CSF_EUR_',
                                 unionSets.dt.uniqID.protFound$Analytes[i], '.regenie.gz'))
  regeniesentOUT.snp <- regeniesentOUT[ID == unionSets.dt.uniqID.protFound$snpID[i]]
  regeniesentOUT.snp$protID.new <- unionSets.dt.uniqID.protFound$Analytes[i]
  regeniesentOUT.snp$protVarBlockID <- unionSets.dt.uniqID.protFound$protVarBlockID[i]
  return(regeniesentOUT.snp)
}


regeniesentOUT.all_T03 # 3696 found --> "Number of ignored tests due to low MAC : 149950"

regeniesentOUT.all_T03.compare <- merge(regeniesentOUT.all_T03, 
                                        unionSets.dt.uniqID.protFound[, list(protVarBlockID, setName, LDname, LP)], 
                                        by = 'protVarBlockID',
                                        all.y = TRUE)
regeniesentOUT.all_T03.compare[is.na(LOG10P)]
regeniesentOUT.all_T03.compare[is.na(LOG10P)]$LOG10P <- 1

regeniesentOUT.all_T03.compare[LDname == '394_START']
cor.test(regeniesentOUT.all_T03.compare[LDname == '394_START']$LP,
         regeniesentOUT.all_T03.compare[LDname == '394_START']$LOG10P)

ggplot(regeniesentOUT.all_T03.compare[LDname == '394_START'], aes(LP, LOG10P)) + geom_point()
ggplot(regeniesentOUT.all_T03.compare[LDname == '394_START'], aes(LP, LOG10P)) + 
  geom_point() +
  xlim(0, 20) +
  ylim(0, 20)



### T04
outdir.T04 <- '/03-DryLab/04-Analyses/2021_Multi-Omics_CC/2022_MultiTissue-MultiOmics_Chengran/t16_benchmarkQTL/run01_CSFlongsLOG10z/d04_regenieOUTnoPEER/t02_sentinelVars/outputS2/'
list.files(outdir.T04)

## to assemble into one table for the downstream comparisons between models
# regeniesentOUT.all_T04 <- foreach(i = 1:10, .combine = 'rbind') %do% {
regeniesentOUT.all_T04 <- foreach(i = 1:nrow(unionSets.dt.uniqID.protFound), .combine = 'rbind') %do% {
  regeniesentOUT <- fread(paste0(outdir.T04, 'T04_WU_CSF_EUR_',
                                unionSets.dt.uniqID.protFound$Analytes[i], '.regenie.gz'))
  regeniesentOUT.snp <- regeniesentOUT[ID == unionSets.dt.uniqID.protFound$snpID[i]]
  regeniesentOUT.snp$protID.new <- unionSets.dt.uniqID.protFound$Analytes[i]
  regeniesentOUT.snp$protVarBlockID <- unionSets.dt.uniqID.protFound$protVarBlockID[i]
  return(regeniesentOUT.snp)
}


regeniesentOUT.all_T04 # 3696 found --> "Number of ignored tests due to low MAC : 149950"

regeniesentOUT.all_T04.compare <- merge(regeniesentOUT.all_T04, 
                                       unionSets.dt.uniqID.protFound[, list(protVarBlockID, setName, LDname, LP)], 
                                       by = 'protVarBlockID',
                                       all.y = TRUE)
regeniesentOUT.all_T04.compare[is.na(LOG10P)]
regeniesentOUT.all_T04.compare[is.na(LOG10P)]$LOG10P <- 1

regeniesentOUT.all_T04.compare[LDname == '394_START']
cor.test(regeniesentOUT.all_T04.compare[LDname == '394_START']$LP,
         regeniesentOUT.all_T04.compare[LDname == '394_START']$LOG10P)

ggplot(regeniesentOUT.all_T04.compare[LDname == '394_START'], aes(LP, LOG10P)) + geom_point()
ggplot(regeniesentOUT.all_T04.compare[LDname == '394_START'], aes(LP, LOG10P)) + 
  geom_point() +
  xlim(0, 20) +
  ylim(0, 20)


### T07
outdir.T07 <- '/03-DryLab/04-Analyses/2021_Multi-Omics_CC/2022_MultiTissue-MultiOmics_Chengran/t16_benchmarkQTL/run02_CSFlongsINVR/d07_regenieOUTwiPEER/t02_sentinelVars/outputS2/'
list.files(outdir.T07)

## to assemble into one table for the downstream comparisons between models
# regeniesentOUT.all_T07 <- foreach(i = 1:10, .combine = 'rbind') %do% {
regeniesentOUT.all_T07 <- foreach(i = 1:nrow(unionSets.dt.uniqID.protFound), .combine = 'rbind') %do% {
  regeniesentOUT <- fread(paste0(outdir.T07, 'T07_WU_CSF_EUR_',
                                 unionSets.dt.uniqID.protFound$Analytes[i], '.regenie.gz'))
  regeniesentOUT.snp <- regeniesentOUT[ID == unionSets.dt.uniqID.protFound$snpID[i]]
  regeniesentOUT.snp$protID.new <- unionSets.dt.uniqID.protFound$Analytes[i]
  regeniesentOUT.snp$protVarBlockID <- unionSets.dt.uniqID.protFound$protVarBlockID[i]
  return(regeniesentOUT.snp)
}


regeniesentOUT.all_T07 # 3696 found --> "Number of ignored tests due to low MAC : 149950"

regeniesentOUT.all_T07.compare <- merge(regeniesentOUT.all_T07, 
                                        unionSets.dt.uniqID.protFound[, list(protVarBlockID, setName, LDname, LP)], 
                                        by = 'protVarBlockID',
                                        all.y = TRUE)
regeniesentOUT.all_T07.compare[is.na(LOG10P)]
regeniesentOUT.all_T07.compare[is.na(LOG10P)]$LOG10P <- 1

regeniesentOUT.all_T07.compare[LDname == '394_START']
cor.test(regeniesentOUT.all_T07.compare[LDname == '394_START']$LP,
         regeniesentOUT.all_T07.compare[LDname == '394_START']$LOG10P)

ggplot(regeniesentOUT.all_T07.compare[LDname == '394_START'], aes(LP, LOG10P)) + geom_point()
ggplot(regeniesentOUT.all_T07.compare[LDname == '394_START'], aes(LP, LOG10P)) + 
  geom_point() +
  xlim(0, 20) +
  ylim(0, 20)

### T08
outdir.T08 <- '/03-DryLab/04-Analyses/2021_Multi-Omics_CC/2022_MultiTissue-MultiOmics_Chengran/t16_benchmarkQTL/run02_CSFlongsINVR/d08_regenieOUTnoPEER/t02_sentinelVars/outputS2/'
list.files(outdir.T08)

## to assemble into one table for the downstream comparisons between models
# regeniesentOUT.all_T08 <- foreach(i = 1:10, .combine = 'rbind') %do% {
regeniesentOUT.all_T08 <- foreach(i = 1:nrow(unionSets.dt.uniqID.protFound), .combine = 'rbind') %do% {
  regeniesentOUT <- fread(paste0(outdir.T08, 'T08_WU_CSF_EUR_',
                                 unionSets.dt.uniqID.protFound$Analytes[i], '.regenie.gz'))
  regeniesentOUT.snp <- regeniesentOUT[ID == unionSets.dt.uniqID.protFound$snpID[i]]
  regeniesentOUT.snp$protID.new <- unionSets.dt.uniqID.protFound$Analytes[i]
  regeniesentOUT.snp$protVarBlockID <- unionSets.dt.uniqID.protFound$protVarBlockID[i]
  return(regeniesentOUT.snp)
}


regeniesentOUT.all_T08 # 3696 found --> "Number of ignored tests due to low MAC : 149950"

regeniesentOUT.all_T08.compare <- merge(regeniesentOUT.all_T08, 
                                        unionSets.dt.uniqID.protFound[, list(protVarBlockID, setName, LDname, LP)], 
                                        by = 'protVarBlockID',
                                        all.y = TRUE)
regeniesentOUT.all_T08.compare[is.na(LOG10P)]
regeniesentOUT.all_T08.compare[is.na(LOG10P)]$LOG10P <- 1

regeniesentOUT.all_T08.compare[LDname == '394_START']
cor.test(regeniesentOUT.all_T08.compare[LDname == '394_START']$LP,
         regeniesentOUT.all_T08.compare[LDname == '394_START']$LOG10P)

ggplot(regeniesentOUT.all_T08.compare[LDname == '394_START'], aes(LP, LOG10P)) + geom_point()
ggplot(regeniesentOUT.all_T08.compare[LDname == '394_START'], aes(LP, LOG10P)) + 
  geom_point() +
  xlim(0, 20) +
  ylim(0, 20)

regeniesentOUT.all_T08.compare[, pqtlTYPE := tstrsplit(protVarBlockID, '_')[4]]
table(regeniesentOUT.all_T08.compare$pqtlTYPE)
regeniesentOUT.all_T08.compare[pqtlTYPE == 'cis']

#################################
### write out the assembled file from each model
new.outdir.T03 <- '/03-DryLab/04-Analyses/2021_Multi-Omics_CC/2022_MultiTissue-MultiOmics_Chengran/t16_benchmarkQTL/run01_CSFlongsLOG10z/d03_regenieOUTwiPEER/t02_sentinelVars/'
list.files(new.outdir.T03)

fwrite(regeniesentOUT.all_T03.compare, paste0(new.outdir.T03, 'assembled_out_T03.csv'))


new.outdir.T04 <- '/03-DryLab/04-Analyses/2021_Multi-Omics_CC/2022_MultiTissue-MultiOmics_Chengran/t16_benchmarkQTL/run01_CSFlongsLOG10z/d04_regenieOUTnoPEER/t02_sentinelVars/'
list.files(new.outdir.T04)

fwrite(regeniesentOUT.all_T04.compare, paste0(new.outdir.T04, 'assembled_out_T04.csv'))

new.outdir.T07 <- '/03-DryLab/04-Analyses/2021_Multi-Omics_CC/2022_MultiTissue-MultiOmics_Chengran/t16_benchmarkQTL/run02_CSFlongsINVR/d07_regenieOUTwiPEER/t02_sentinelVars/'
list.files(new.outdir.T07)

fwrite(regeniesentOUT.all_T07.compare, paste0(new.outdir.T07, 'assembled_out_T07.csv'))

new.outdir.T08 <- '/03-DryLab/04-Analyses/2021_Multi-Omics_CC/2022_MultiTissue-MultiOmics_Chengran/t16_benchmarkQTL/run02_CSFlongsINVR/d08_regenieOUTnoPEER/t02_sentinelVars/'
list.files(new.outdir.T08)

fwrite(regeniesentOUT.all_T08.compare, paste0(new.outdir.T08, 'assembled_out_T08.csv'))























