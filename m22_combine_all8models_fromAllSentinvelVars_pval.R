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

#################################
### load the assembled files from four models from plink2
new.outdir.T01 <- '/03-DryLab/04-Analyses/2021_Multi-Omics_CC/2022_MultiTissue-MultiOmics_Chengran/t16_benchmarkQTL/run01_CSFlongsLOG10z/d01_plink2OUTwiPEER/t02_sentinelVars/'
list.files(new.outdir.T01)
model.out.T01 <- fread(paste0(new.outdir.T01, 'assembled_out_T01.csv'))

new.outdir.T02 <- '/03-DryLab/04-Analyses/2021_Multi-Omics_CC/2022_MultiTissue-MultiOmics_Chengran/t16_benchmarkQTL/run01_CSFlongsLOG10z/d02_plink2OUTnoPEER/t02_sentinelVars/'
list.files(new.outdir.T02)
model.out.T02 <- fread(paste0(new.outdir.T02, 'assembled_out_T02.csv'))

new.outdir.T05 <- '/03-DryLab/04-Analyses/2021_Multi-Omics_CC/2022_MultiTissue-MultiOmics_Chengran/t16_benchmarkQTL/run02_CSFlongsINVR/d05_plink2OUTwiPEER/t02_sentinelVars/'
list.files(new.outdir.T05)
model.out.T05 <- fread(paste0(new.outdir.T05, 'assembled_out_T05.csv'))

new.outdir.T06 <- '/03-DryLab/04-Analyses/2021_Multi-Omics_CC/2022_MultiTissue-MultiOmics_Chengran/t16_benchmarkQTL/run02_CSFlongsINVR/d06_plink2OUTnoPEER/t02_sentinelVars/'
list.files(new.outdir.T06)
model.out.T06 <- fread(paste0(new.outdir.T06, 'assembled_out_T06.csv'))

### combine all four tables into one
all(model.out.T01$protVarBlockID == model.out.T02$protVarBlockID)
all(model.out.T02$protVarBlockID == model.out.T05$protVarBlockID)
all(model.out.T05$protVarBlockID == model.out.T06$protVarBlockID)


model1n2.combine.plink2.dt <- merge(model.out.T01[, list(protVarBlockID, LOG10_P)],
                                    model.out.T02[, list(protVarBlockID, LOG10_P)],
                                    by = 'protVarBlockID', suffixes = c('.T01', '.T02'))
model5n6.combine.plink2.dt <- merge(model.out.T05[, list(protVarBlockID, LOG10_P)],
                                    model.out.T06[, list(protVarBlockID, LOG10_P)],
                                    by = 'protVarBlockID', suffixes = c('.T05', '.T06'))
modelAll4.combine.plink2.dt <- merge(model1n2.combine.plink2.dt, model5n6.combine.plink2.dt, by = 'protVarBlockID')
modelAll4.combine.plink2.dt


summary(as.matrix(modelAll4.combine.plink2.dt[, 2:5]))
modelAll4.combine.plink2.dt[is.na(LOG10_P.T01)]$LOG10_P.T01 <- 0
modelAll4.combine.plink2.dt[is.na(LOG10_P.T02)]$LOG10_P.T02 <- 0
modelAll4.combine.plink2.dt[is.na(LOG10_P.T05)]$LOG10_P.T05 <- 0
modelAll4.combine.plink2.dt[is.na(LOG10_P.T06)]$LOG10_P.T06 <- 0
# library(ComplexHeatmap)
# Heatmap(as.matrix(modelAll4.combine.plink2.dt[1:100, 2:5]))

# long.modelAll4.combine.plink2.dt <- data.table::melt(modelAll4.combine.plink2.dt, id.vars = 1)
# long.modelAll4.combine.plink2.dt
# long.modelAll4.combine.plink2.dt[, loglogP := log10(value)]
# long.modelAll4.combine.plink2.dt[is.na(loglogP)]
# ggplot(long.modelAll4.combine.plink2.dt[protVarBlockID %in% 
#                                           modelAll4.combine.plink2.dt[1:100]$protVarBlockID], 
#        aes(variable, protVarBlockID)) +
#   geom_raster(aes(fill = loglogP))
# 
# ggplot(long.modelAll4.combine.plink2.dt[protVarBlockID %in% 
#                                           modelAll4.combine.plink2.dt[1:100]$protVarBlockID], 
#        aes(variable, protVarBlockID)) +
#   geom_raster(aes(fill = value))

#################################
### load the assembled files from four models from regenie
new.outdir.T03 <- '/03-DryLab/04-Analyses/2021_Multi-Omics_CC/2022_MultiTissue-MultiOmics_Chengran/t16_benchmarkQTL/run01_CSFlongsLOG10z/d03_regenieOUTwiPEER/t02_sentinelVars/'
list.files(new.outdir.T03)
model.out.T03 <- fread(paste0(new.outdir.T03, 'assembled_out_T03.csv'))

new.outdir.T04 <- '/03-DryLab/04-Analyses/2021_Multi-Omics_CC/2022_MultiTissue-MultiOmics_Chengran/t16_benchmarkQTL/run01_CSFlongsLOG10z/d04_regenieOUTnoPEER/t02_sentinelVars/'
list.files(new.outdir.T04)
model.out.T04 <- fread(paste0(new.outdir.T04, 'assembled_out_T04.csv'))

new.outdir.T07 <- '/03-DryLab/04-Analyses/2021_Multi-Omics_CC/2022_MultiTissue-MultiOmics_Chengran/t16_benchmarkQTL/run02_CSFlongsINVR/d07_regenieOUTwiPEER/t02_sentinelVars/'
list.files(new.outdir.T07)
model.out.T07 <- fread(paste0(new.outdir.T07, 'assembled_out_T07.csv'))

new.outdir.T08 <- '/03-DryLab/04-Analyses/2021_Multi-Omics_CC/2022_MultiTissue-MultiOmics_Chengran/t16_benchmarkQTL/run02_CSFlongsINVR/d08_regenieOUTnoPEER/t02_sentinelVars/'
list.files(new.outdir.T08)
model.out.T08 <- fread(paste0(new.outdir.T08, 'assembled_out_T08.csv'))

### combine all four tables into one
all(model.out.T03$protVarBlockID == model.out.T04$protVarBlockID)
all(model.out.T04$protVarBlockID == model.out.T07$protVarBlockID)
all(model.out.T07$protVarBlockID == model.out.T08$protVarBlockID)


model3n4.combine.regenie.dt <- merge(model.out.T03[, list(protVarBlockID, LOG10_P=LOG10P)],
                                    model.out.T04[, list(protVarBlockID, LOG10_P=LOG10P)],
                                    by = 'protVarBlockID', suffixes = c('.T03', '.T04'))
model7n8.combine.regenie.dt <- merge(model.out.T07[, list(protVarBlockID, LOG10_P=LOG10P)],
                                    model.out.T08[, list(protVarBlockID, LOG10_P=LOG10P)],
                                    by = 'protVarBlockID', suffixes = c('.T07', '.T08'))
modelAll4.combine.regenie.dt <- merge(model3n4.combine.regenie.dt, model7n8.combine.regenie.dt, by = 'protVarBlockID')
modelAll4.combine.regenie.dt


summary(as.matrix(modelAll4.combine.regenie.dt[, 2:5]))
# modelAll4.combine.regenie.dt[is.na(LOG10_P.T03)]$LOG10_P.T03 <- 0
# modelAll4.combine.regenie.dt[is.na(LOG10_P.T04)]$LOG10_P.T04 <- 0
# modelAll4.combine.regenie.dt[is.na(LOG10_P.T07)]$LOG10_P.T07 <- 0
# modelAll4.combine.regenie.dt[is.na(LOG10_P.T08)]$LOG10_P.T08 <- 0

#################################
### combined plink2 and regenie for all 8 models
modelAll8.dt <- merge(modelAll4.combine.plink2.dt, 
                      modelAll4.combine.regenie.dt, by = 'protVarBlockID')


long.modelAll8.dt <- data.table::melt(modelAll8.dt, id.vars = 1)

#################################
#### check the genome-wide hits per model and rank across 8 models
long.modelAll8.dt %>%
  filter(value > 7.3) %>%
  group_by(variable) %>%
  summarise(count_pqtl = n()) %>%
  arrange(-count_pqtl)

## as expected, all 4 peer+ models rank higher than all 4 peer- models  
## also, plink2 models higher than regenie with the same-covar+same-prot,e.g. T05 > T07; T01 > T03
## also, INVR models higher than log10z with the same-tool+same-prot,e.g. T05 > T01; T07 > T03

## draw the upsetR plot
library(UpSetR)
upset(fromList(list(T01=long.modelAll8.dt[value > 7.3 & variable == 'LOG10_P.T01']$protVarBlockID,
                    T02=long.modelAll8.dt[value > 7.3 & variable == 'LOG10_P.T02']$protVarBlockID,
                    T03=long.modelAll8.dt[value > 7.3 & variable == 'LOG10_P.T03']$protVarBlockID,
                    T04=long.modelAll8.dt[value > 7.3 & variable == 'LOG10_P.T04']$protVarBlockID,
                    T05=long.modelAll8.dt[value > 7.3 & variable == 'LOG10_P.T05']$protVarBlockID,
                    T06=long.modelAll8.dt[value > 7.3 & variable == 'LOG10_P.T06']$protVarBlockID,
                    T07=long.modelAll8.dt[value > 7.3 & variable == 'LOG10_P.T07']$protVarBlockID,
                    T08=long.modelAll8.dt[value > 7.3 & variable == 'LOG10_P.T08']$protVarBlockID)),
      order.by = "freq", nsets = 8)





### split by cis and trans qtlTYPE
long.modelAll8.dt[, qtlTYPE := tstrsplit(protVarBlockID, '_')[4]]
long.modelAll8.dt
# long.modelAll8.dt %>%
#   filter(value > 7.3) %>%
#   group_by(variable, qtlTYPE) %>%
#   summarise(count_pqtl = n()) %>%
#   arrange(-count_pqtl)


long.modelAll8.dt[qtlTYPE == 'cis'] %>%
  filter(value > 7.3) %>%
  group_by(variable) %>%
  summarise(count_pqtl = n()) %>%
  arrange(-count_pqtl)

long.modelAll8.dt[qtlTYPE == 'trans'] %>%
  filter(value > 7.3) %>%
  group_by(variable) %>%
  summarise(count_pqtl = n()) %>%
  arrange(-count_pqtl)

long.modelAll8.dt[qtlTYPE == 'trans'] %>%
  filter(value > 11.1) %>%
  group_by(variable) %>%
  summarise(count_pqtl = n()) %>%
  arrange(-count_pqtl)


## draw the upsetR plot
upset(fromList(list(T01=long.modelAll8.dt[qtlTYPE == 'cis' & value > 7.3 & variable == 'LOG10_P.T01']$protVarBlockID,
                    T02=long.modelAll8.dt[qtlTYPE == 'cis' & value > 7.3 & variable == 'LOG10_P.T02']$protVarBlockID,
                    T03=long.modelAll8.dt[qtlTYPE == 'cis' & value > 7.3 & variable == 'LOG10_P.T03']$protVarBlockID,
                    T04=long.modelAll8.dt[qtlTYPE == 'cis' & value > 7.3 & variable == 'LOG10_P.T04']$protVarBlockID,
                    T05=long.modelAll8.dt[qtlTYPE == 'cis' & value > 7.3 & variable == 'LOG10_P.T05']$protVarBlockID,
                    T06=long.modelAll8.dt[qtlTYPE == 'cis' & value > 7.3 & variable == 'LOG10_P.T06']$protVarBlockID,
                    T07=long.modelAll8.dt[qtlTYPE == 'cis' & value > 7.3 & variable == 'LOG10_P.T07']$protVarBlockID,
                    T08=long.modelAll8.dt[qtlTYPE == 'cis' & value > 7.3 & variable == 'LOG10_P.T08']$protVarBlockID)),
      order.by = "freq", nsets = 8)

upset(fromList(list(T01=long.modelAll8.dt[qtlTYPE == 'trans' & value > 7.3 & variable == 'LOG10_P.T01']$protVarBlockID,
                    T02=long.modelAll8.dt[qtlTYPE == 'trans' & value > 7.3 & variable == 'LOG10_P.T02']$protVarBlockID,
                    T03=long.modelAll8.dt[qtlTYPE == 'trans' & value > 7.3 & variable == 'LOG10_P.T03']$protVarBlockID,
                    T04=long.modelAll8.dt[qtlTYPE == 'trans' & value > 7.3 & variable == 'LOG10_P.T04']$protVarBlockID,
                    T05=long.modelAll8.dt[qtlTYPE == 'trans' & value > 7.3 & variable == 'LOG10_P.T05']$protVarBlockID,
                    T06=long.modelAll8.dt[qtlTYPE == 'trans' & value > 7.3 & variable == 'LOG10_P.T06']$protVarBlockID,
                    T07=long.modelAll8.dt[qtlTYPE == 'trans' & value > 7.3 & variable == 'LOG10_P.T07']$protVarBlockID,
                    T08=long.modelAll8.dt[qtlTYPE == 'trans' & value > 7.3 & variable == 'LOG10_P.T08']$protVarBlockID)),
      order.by = "freq", nsets = 8)
#################################
#### check the genome-wide hits per model and rank across 8 models only on the chr3_LD394-locus
long.modelAll8.dt[grep('394_START_trans', protVarBlockID)] %>%
  filter(value > 7.3) %>%
  group_by(variable) %>%
  summarise(count_pqtl = n()) %>%
  arrange(-count_pqtl)

# long.modelAll8.dt[grep('394_START_trans', protVarBlockID)] %>%
#   filter(value > 11) %>%
#   group_by(variable) %>%
#   summarise(count_pqtl = n()) %>%
#   arrange(-count_pqtl)

## as expected, the transProteins associated with chr3_LD394 only showed up in the peer- models
## T02 is the strongest, as per DW2022 model
## for study-wide significance, maybe the sample size now dropped PPMI, leading to p-value decreasing


#################################
### sort the p-value and check which models has the largest significant findings
setorder(long.modelAll8.dt, -value)

long.modelAll8.dt.uniqID <- unique(long.modelAll8.dt, by = 'protVarBlockID')
nrow(long.modelAll8.dt.uniqID) # 3748

table(long.modelAll8.dt.uniqID$variable)
##T05 has the largest set, but this can be due to the inflation of plink2 with PEER we observed earlier for TREM2 findings



### split trans into cis-found and cis-missing given the same protID
long.modelAll8.dt[, protID := tstrsplit(protVarBlockID, '_')[1]]
long.modelAll8.dt
long.modelAll8.dt.cis <- long.modelAll8.dt[qtlTYPE == 'cis' & value > 7.3]
long.modelAll8.dt.trans <- long.modelAll8.dt[qtlTYPE == 'trans' & value > 7.3]
# long.modelAll8.dt.trans <- long.modelAll8.dt[qtlTYPE == 'trans' & value > 11.1]

table(long.modelAll8.dt.cis$variable)
table(long.modelAll8.dt.trans$variable)

long.modelAll8.dt.trans.cisFound <- long.modelAll8.dt.trans[protID %in% long.modelAll8.dt.cis$protID]
long.modelAll8.dt.trans.cisMissing <- long.modelAll8.dt.trans[!protID %in% long.modelAll8.dt.cis$protID]

table(long.modelAll8.dt.trans.cisFound$variable)
table(long.modelAll8.dt.trans.cisMissing$variable)

## I should use c('variable', 'protID') rather than protID to check the uniqueness of proteins!!
long.modelAll8.dt.trans.cisFound.uniq <- unique(long.modelAll8.dt.trans.cisFound, by = c('variable', 'protID'))
long.modelAll8.dt.trans.cisMissing.uniq <- unique(long.modelAll8.dt.trans.cisMissing, by = c('variable', 'protID'))

long.modelAll8.dt.trans.cisFound.uniq %>%
  group_by(variable) %>%
  summarise(count_pqtl = n()) %>%
  arrange(-count_pqtl) -> check_trans_cisFound.tb

long.modelAll8.dt.trans.cisMissing.uniq %>%
  group_by(variable) %>%
  summarise(count_pqtl = n()) %>%
  arrange(-count_pqtl) -> check_trans_cisMissing.tb


trans.gw.compare <- data.table(merge(check_trans_cisFound.tb, check_trans_cisMissing.tb, 
      by = 'variable', suffixes = c('.cisFound', '.cisMissing')))

trans.gw.compare
trans.gw.compare[, trans_sum := count_pqtl.cisFound+count_pqtl.cisMissing]
trans.gw.compare[, ratio.Found := count_pqtl.cisFound/trans_sum]
setorder(trans.gw.compare, -ratio.Found)
trans.gw.compare # T07 ranked higher than T05, and all other 7 models

trans.gw.compare[, ratio.Odds := count_pqtl.cisFound/count_pqtl.cisMissing]
setorder(trans.gw.compare, -ratio.Odds)
trans.gw.compare
