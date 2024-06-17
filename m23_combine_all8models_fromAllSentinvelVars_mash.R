rm(list = ls())
library(dplyr)
library(data.table)
library(reshape2)
library(openxlsx)
library(Biobase)
library(ggplot2)
library(cowplot)

library(foreach)
library(mashr)
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


model1n2.combine.plink2.dt <- merge(model.out.T01[, list(protVarBlockID, BETA, SE)],
                                    model.out.T02[, list(protVarBlockID, BETA, SE)],
                                    by = 'protVarBlockID', suffixes = c('.T01', '.T02'))
model5n6.combine.plink2.dt <- merge(model.out.T05[, list(protVarBlockID, BETA, SE)],
                                    model.out.T06[, list(protVarBlockID, BETA, SE)],
                                    by = 'protVarBlockID', suffixes = c('.T05', '.T06'))
modelAll4.combine.plink2.dt <- merge(model1n2.combine.plink2.dt, model5n6.combine.plink2.dt, by = 'protVarBlockID')
modelAll4.combine.plink2.dt


summary(as.matrix(modelAll4.combine.plink2.dt[, 2:5]))
modelAll4.combine.plink2.dt[is.na(BETA.T01)]$BETA.T01 <- 0
modelAll4.combine.plink2.dt[is.na(BETA.T02)]$BETA.T02 <- 0
modelAll4.combine.plink2.dt[is.na(BETA.T05)]$BETA.T05 <- 0
modelAll4.combine.plink2.dt[is.na(BETA.T06)]$BETA.T06 <- 0

modelAll4.combine.plink2.dt[is.na(SE.T01)]$SE.T01 <- 1e6
modelAll4.combine.plink2.dt[is.na(SE.T02)]$SE.T02 <- 1e6
modelAll4.combine.plink2.dt[is.na(SE.T05)]$SE.T05 <- 1e6
modelAll4.combine.plink2.dt[is.na(SE.T06)]$SE.T06 <- 1e6


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


model3n4.combine.regenie.dt <- merge(model.out.T03[, list(protVarBlockID, BETA, SE)],
                                    model.out.T04[, list(protVarBlockID, BETA, SE)],
                                    by = 'protVarBlockID', suffixes = c('.T03', '.T04'))
model7n8.combine.regenie.dt <- merge(model.out.T07[, list(protVarBlockID, BETA, SE)],
                                    model.out.T08[, list(protVarBlockID, BETA, SE)],
                                    by = 'protVarBlockID', suffixes = c('.T07', '.T08'))
modelAll4.combine.regenie.dt <- merge(model3n4.combine.regenie.dt, model7n8.combine.regenie.dt, by = 'protVarBlockID')
modelAll4.combine.regenie.dt


summary(as.matrix(modelAll4.combine.regenie.dt[, 2:5]))
modelAll4.combine.regenie.dt[is.na(BETA.T03)]$BETA.T03 <- 0
modelAll4.combine.regenie.dt[is.na(BETA.T04)]$BETA.T04 <- 0
modelAll4.combine.regenie.dt[is.na(BETA.T07)]$BETA.T07 <- 0
modelAll4.combine.regenie.dt[is.na(BETA.T08)]$BETA.T08 <- 0

modelAll4.combine.regenie.dt[is.na(SE.T03)]$SE.T03 <- 1e6
modelAll4.combine.regenie.dt[is.na(SE.T04)]$SE.T04 <- 1e6
modelAll4.combine.regenie.dt[is.na(SE.T07)]$SE.T07 <- 1e6
modelAll4.combine.regenie.dt[is.na(SE.T08)]$SE.T08 <- 1e6


#################################
### combined plink2 and regenie for all 8 models
modelAll8.dt <- merge(modelAll4.combine.plink2.dt, 
                      modelAll4.combine.regenie.dt, by = 'protVarBlockID')


#################################
## perform mash analysis
#### step-1: Read in the data
input.data02 = mash_set_data(as.matrix(modelAll8.dt[, list(BETA.T01, BETA.T02, BETA.T03, BETA.T04,
                                                           BETA.T05, BETA.T06, BETA.T07, BETA.T08)]), 
                             as.matrix(modelAll8.dt[, list(SE.T01, SE.T02, SE.T03, SE.T04,
                                                           SE.T05, SE.T06, SE.T07, SE.T08)]))
#### step-2: Set up the covariance matrices
U.c02 = cov_canonical(input.data02)  
print(names(U.c02))

#### Step 3: fit the model
m.c02 = mash(input.data02, U.c02)

#### Step 4: Extract Posterior Summaries
head(get_lfsr(m.c02))

lfsr.mat02 <- get_lfsr(m.c02)
lfsr.dt02 <- data.table(lfsr.mat02)

lfsr.dt02[, protVarBlockID := modelAll8.dt$protVarBlockID]


##### sharing
print(get_pairwise_sharing(m.c02))
summary(get_pairwise_sharing(m.c02, factor = 0.5)) # 70% sharing given default factor as 0.5
library(ComplexHeatmap)
Heatmap(get_pairwise_sharing(m.c02, factor = 0.5))



print(get_pairwise_sharing(m.c02, factor = 0.1))
summary(get_pairwise_sharing(m.c02, factor = 0.1)) # 73% sharing given factor as 0.1
print(get_pairwise_sharing(m.c02, factor=0))

# print(get_pairwise_sharing(m.c02, factor = 0.1, FUN = abs))
#### this is to dissect how sharing % is coming from
postMean.dt02 <- data.table(m.c02$result$PosteriorMean)
postMean.dt02

# postMean.dt02[, factor := postMean_AFR/postMean_EUR]
# postMean.dt02
# postMean.dt02[, pQTL := compare_effAFR.dt$pQTL]
# postMean.dt02[, share := (factor > 1/10 & factor < 10)]
# table(postMean.dt02$share)
# table(postMean.dt02$share)/nrow(postMean.dt02)



