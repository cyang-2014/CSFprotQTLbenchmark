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


model1n2.combine.plink2.dt <- merge(model.out.T01[, list(protVarBlockID, chr.var=`#CHROM`, pos.var=POS, LOG10_P)],
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
modelAll8.dt

long.modelAll8.dt <- data.table::melt(modelAll8.dt, id.vars = 1:3)

#################################
#### check the genome-wide hits per model and rank across 8 models
long.modelAll8.dt %>%
  filter(value > 7.3) %>%
  group_by(variable) %>%
  summarise(count_pqtl = n()) %>%
  arrange(-count_pqtl)





###########################################
## to draw the 2d manhattan plot, I need to load in gene-position file for CSF-prot-data, not the plasma-prot-data
DIRsoma7kgeneLoc <- '/03-DryLab/04-Analyses/2021_Multi-Omics_CC/2022_MultiTissue-MultiOmics_Chengran/t07_CSFprotMay2023/run01_1Longs2Ppmi3Stanford_rmFACE/d03_input5June2023/'
list.files(DIRsoma7kgeneLoc)

LONGSsoma7k.feature.geneCoord <- fread(paste0(DIRsoma7kgeneLoc, 'f06_allSplitFeature7390human_geneCoordinates.csv'))
LONGSsoma7k.feature.geneCoord

LONGSsoma7k.feature.geneCoord[, Analytes := tstrsplit(AnalytesUniprotSplit, '_')[1]]

### change the "chr" into integer class instead of character, for the dot-plot drawing
genepos.input.chrAUTO <- LONGSsoma7k.feature.geneCoord[!chr %in% c('chrX', 'chrY')]
genepos.input.chrX <- LONGSsoma7k.feature.geneCoord[chr %in% c('chrX')]
genepos.input.chrY <- LONGSsoma7k.feature.geneCoord[chr %in% c('chrY')]


genepos.input.chrAUTO[, chr := tstrsplit(chr, 'chr')[2]]
genepos.input.chrAUTO[, chr := as.integer(chr)]

genepos.input.chrX[, chr := 23]
genepos.input.chrY[, chr := 24]

genepos.input.newCHR <- rbindlist(list(genepos.input.chrAUTO, genepos.input.chrX, genepos.input.chrY))
nrow(genepos.input.newCHR)
genepos.input.newCHR[, chr := as.integer(chr)]

# genepos.input.newCHR.uniq <- unique(genepos.input.newCHR, by = 'Analytes', fromLast = FALSE)
# nrow(genepos.input.newCHR.uniq)

#####################################
## add protein-encoded gene-position info
modelAll8.dt[, Analytes := tstrsplit(protVarBlockID, '_')[1]]
modelAll8.dt
modelAll8.dt[, type.pre := tstrsplit(protVarBlockID, '_')[4]]
table(modelAll8.dt$type.pre)


modelAll8.dt.withGenePos <- data.table(left_join(modelAll8.dt,
                                                 genepos.input.newCHR[, list(Analytes, 
                                                                                   chr.prot=chr, 
                                                                                   start.prot=start)],
                                                 by = 'Analytes'))
modelAll8.dt.withGenePos


modelAll8.dt.withGenePos[type.pre == 'cis'][!chr.var == chr.prot]
modelAll8.dt.withGenePos[Analytes == 'X10365.132']

modelAll8.dt.withGenePos.cis <- modelAll8.dt.withGenePos[type.pre == 'cis']
modelAll8.dt.withGenePos.trans <- modelAll8.dt.withGenePos[type.pre == 'trans']

## keep the unique pairs after determining cis and trans-pQTL by the SNP and gene position
modelAll8.dt.withGenePos.cis.uniq <- modelAll8.dt.withGenePos.cis[chr.var==chr.prot]
modelAll8.dt.withGenePos.cis.uniq





modelAll8.dt.withGenePos.trans[, type.check := ifelse(chr.var != chr.prot, 'trans',
                                                         ifelse((pos.var > start.prot + 1e6) |
                                                                  (pos.var < start.prot - 1e6),
                                                                'trans', 'cis'))]
modelAll8.dt.withGenePos.trans[type.check == 'cis']
modelAll8.dt.withGenePos.trans.uniq <- modelAll8.dt.withGenePos.trans[!type.check == 'cis']
modelAll8.dt.withGenePos.trans.uniq[, type.check := NULL]

modelAll8.dt.withGenePos.correct <- rbindlist(list(modelAll8.dt.withGenePos.cis.uniq, modelAll8.dt.withGenePos.trans.uniq))

modelAll8.dt.withGenePos.correct[duplicated(protVarBlockID)]
modelAll8.dt.withGenePos.correct[protVarBlockID == 'X18821.9_MHC_START_cis']

modelAll8.dt.withGenePos.correct.uniq <- unique(modelAll8.dt.withGenePos.correct, by = 'protVarBlockID')
nrow(modelAll8.dt.withGenePos.correct.uniq)
table(modelAll8.dt.withGenePos.correct.uniq$type.pre)

##################################################################
## draw the dot plot to visualize the pQTLs in one figure

# format df for SNP
df.tmp.SNP <- modelAll8.dt.withGenePos.correct.uniq %>% 
  
  # Compute chromosome size
  group_by(chr.var) %>% 
  summarise(chr_len=max(pos.var)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot_SNP =cumsum(as.numeric(chr_len))-chr_len) 

# %>%
#   select(-chr_len) 


# Add this info to the initial dataset
df.tmp2.SNP <- left_join(modelAll8.dt.withGenePos.correct.uniq, df.tmp.SNP, by=c("chr.var")) 

# Add a cumulative position of each SNP
df.tmp3.SNP <- df.tmp2.SNP %>%
  arrange(chr.var, pos.var) %>%
  mutate(BPcum_SNP=pos.var+tot_SNP)

# get chromosome center positions for x-axis
axisdf_SNP <- df.tmp3.SNP %>% group_by(chr.var) %>% summarize(center=( max(BPcum_SNP) + min(BPcum_SNP) ) / 2,
                                                            boundary = max(BPcum_SNP))

# format df for protein
df.tmp.protein <- modelAll8.dt.withGenePos %>% 
  
  # Compute chromosome size
  group_by(chr.prot) %>% 
  summarise(chr_len=max(start.prot)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot_protein =cumsum(as.numeric(chr_len))-chr_len) 
# %>%
#   select(-chr_len) 


# Add this info to the initial dataset
df.tmp2.protein <- left_join(df.tmp3.SNP, df.tmp.protein, by=c("chr.prot")) 

# Add a cumulative position of each protein
df.tmp3.protein <- df.tmp2.protein %>%
  arrange(chr.prot, start.prot) %>%
  mutate( BPcum_protein=start.prot+tot_protein)

# get chromosome center positions for x-axis
axisdf_protein <- df.tmp3.protein %>% group_by(chr.prot) %>% summarize(center=( max(BPcum_protein) + min(BPcum_protein) ) / 2,
                                                                    boundary = max(BPcum_protein))



## generate the 2d-manhattan (dot-plot)
df.tmp3.protein
df.tmp3.protein[, type.pqtl := tstrsplit(protVarBlockID, '_')[4]]
### now the cis pqtl issue with "inter-chrom" info was resolved on 4/1/2024
df.tmp3.protein[type.pqtl == 'cis'][!chr.var == chr.prot]
# genepos.input.newCHR[Analytes == 'X16927.9']
table(df.tmp3.protein$type.pqtl, df.tmp3.protein$type.pre)



# ### color code by cis and trans
# ggplot(df.tmp3.protein, aes(x=BPcum_SNP, y=BPcum_protein)) +
#   # Show all points
#   geom_point(aes(color=as.factor(type.pqtl)), alpha=0.8, size=1) +
#   scale_color_manual(values = c('magenta', 'dodgerblue')) +
#   # custom X axis:
#   scale_x_continuous(label = axisdf_SNP$chr.var, 
#                      breaks= axisdf_SNP$center,
#                      minor_breaks = axisdf_SNP$boundary) +
#   scale_y_continuous(label = axisdf_protein$chr.prot, 
#                      breaks= axisdf_protein$center,
#                      minor_breaks = axisdf_protein$boundary)  +
#   # add plot and axis titles
#   ggtitle(paste0("Genomic locations of pQTLs")) +
#   labs(x = "pQTL position", y = 'Protein position') +
#   theme_cowplot(font_size = 8) +
#   theme(axis.text.x = element_text(angle = 60), 
#         legend.position = "bottom",
#         panel.grid.minor.x = element_line(colour="black", size=0.5),
#         panel.grid.minor.y = element_line(colour="black", size=0.5))


##################################################################################
### color code by model-ID after melting into long form and subset per p-value for GW hits within each model only
df.tmp3.protein.long <- data.table::melt(df.tmp3.protein[, list(protVarBlockID, BPcum_SNP, BPcum_protein, type.pqtl,
                                                    LOG10_P.T01, LOG10_P.T02, LOG10_P.T03, LOG10_P.T04,
                                                    LOG10_P.T05, LOG10_P.T06, LOG10_P.T07, LOG10_P.T08)], 
                                         id.vars = 1:4, variable.name = 'modelID', value.name = 'LOG10_P')
df.tmp3.protein.long
  
ggplot(df.tmp3.protein.long[LOG10_P > -log10(5e-8)], aes(x=BPcum_SNP, y=BPcum_protein)) +
  # Show all points
  geom_point(aes(color=as.factor(modelID)), alpha=0.8, size=1) +
  scale_color_manual(values=c("darkred", "red", "orange", "yellow",
                              "lightgreen", "green", "lightblue", "blue"))+
  # custom X axis:
  scale_x_continuous(label = axisdf_SNP$chr.var, 
                     breaks= axisdf_SNP$center,
                     minor_breaks = axisdf_SNP$boundary) +
  scale_y_continuous(label = axisdf_protein$chr.prot, 
                     breaks= axisdf_protein$center,
                     minor_breaks = axisdf_protein$boundary)  +
  # add plot and axis titles
  ggtitle(paste0("Genomic locations of pQTLs")) +
  labs(x = "pQTL position", y = 'Protein position') +
  theme_cowplot(font_size = 8) +
  theme(axis.text.x = element_text(angle = 60), 
        legend.position = "bottom",
        panel.grid.minor.x = element_line(colour="gray", size=0.5),
        panel.grid.minor.y = element_line(colour="gray", size=0.5)) +
  facet_wrap(~modelID, ncol = 2)


##################################################################################
## only include two models per Carlos's feedback on 5/6/2024 as the panels for Figure-1
## also in one column rather than 2
## per Ciyang's suggestion, I picked two color-contrast models as T07 and T04
## change gray further to gray95?
df.tmp3.protein.long[LOG10_P > -log10(5e-8) & modelID %in% c('LOG10_P.T07', 'LOG10_P.T04')] %>%
  ggplot(aes(x=BPcum_SNP, y=BPcum_protein)) +
    # Show all points
    geom_point(aes(color=as.factor(modelID)), alpha=0.8, size=1) +
    scale_color_manual(values=c("yellow", "lightblue"))+
    # custom X axis:
    scale_x_continuous(label = axisdf_SNP$chr.var, 
                       breaks= axisdf_SNP$center,
                       minor_breaks = axisdf_SNP$boundary) +
    scale_y_continuous(label = axisdf_protein$chr.prot, 
                       breaks= axisdf_protein$center,
                       minor_breaks = axisdf_protein$boundary)  +
    # add plot and axis titles
    ggtitle(paste0("Genomic locations of pQTLs")) +
    labs(x = "pQTL position", y = 'Protein position') +
    theme_cowplot(font_size = 8) +
    theme(axis.text.x = element_text(angle = 60), 
          legend.position = "bottom",
          panel.grid.minor.x = element_line(colour="gray95", size=0.5),
          panel.grid.minor.y = element_line(colour="gray95", size=0.5)) +
    facet_wrap(~modelID, ncol = 1)

ggsave(filename = 'figs/GWAS_dotplot/CSFbenchmark8models_model04and07only.pdf', width = 10, height = 6)



