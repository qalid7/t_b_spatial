#libraries 
#############.
library(tidyr)
library(plyr)
library(dplyr)
library(reshape2)
library(ggpubr)
library(ggplot2)
library("survival")
library("survminer")
library(tidyverse)
library(RTCGA)
library(RTCGA.clinical)
library(gridExtra)
library(RColorBrewer)
library(grid)
library(ggcorrplot)
library(EDASeq)
library(TCGAbiolinks)
library(pheatmap)
#############.
setwd('/Users/hzhang/Documents/report/lusc_b/scripts/new')

#Functions ####
#############.
plt_hs_celltype <- function(df, compare_hs, cell_type = 'cd8', cell_name = 'CD8+ T cells'){
  p <- ggplot(subset(df, variable == cell_type), aes(x= Hotspots, y = value)) +
    geom_boxplot(outlier.shape = NA, aes(fill = Hotspots))+
    geom_dotplot(binaxis='y', stackdir='center', stackratio = 0.8,
                 position=position_jitter(height = 2, width = 0.1, seed = 233),
                 dotsize=0.7, aes(fill = Hotspots), color = 'black',stroke=0.8) + 
    scale_fill_manual(values= c("#E84D3F","#77BFD1","#DEA15B"))+
    scale_color_manual(values= c("#E84D3F","#77BFD1","#DEA15B"))+
    ylim(c(0,100)) + ylab(paste0('% of ', cell_name)) +
    theme_bw() + 
    theme(axis.text.x = element_blank(), legend.position = 'None', axis.title=element_text(size=18),
          axis.text=element_text(size=16, colour="black"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.background = element_blank(), aspect.ratio = 1.6/1.3,
          axis.ticks.length=unit(.25, "cm")) + xlab('') +
    annotate(geom = 'segment', x = 1, xend = 2, y = 60, yend = 60, color = subset(compare_hs, variable == cell_type)$line_color[1])+
    geom_text(x = 1.5, y = 63, label = subset(compare_hs, variable == cell_type)$p_adj_star[1]) +
    annotate(geom = 'segment', x = 1, xend = 3, y = 75, yend = 75, color = subset(compare_hs, variable == cell_type)$line_color[2])+
    geom_text(x = 2, y = 78, label = subset(compare_hs, variable == cell_type)$p_adj_star[2]) +
    annotate(geom = 'segment', x = 2, xend = 3, y = 90, yend = 90, color = subset(compare_hs, variable == cell_type)$line_color[3])+
    geom_text(x = 2.5, y = 93, label = subset(compare_hs, variable == cell_type)$p_adj_star[3])
  return(p)
}

plt_region_celltype <- function(df, compare_region, cell_type = 'cd8', cell_name = 'CD8+ T cells'){
  p <- ggplot(subset(df_m,Hotspots %in% c('Immune') & variable == cell_type), aes(x= Region, y = value)) +
    geom_boxplot(outlier.shape = NA, aes(fill = Region))+
    geom_dotplot(binaxis='y', stackdir='center', stackratio = 0.3,
                 position=position_jitter(height = 2, width = 0.1, seed = 233),
                 dotsize=0.7, aes(fill = Region), color = 'black',stroke=0.8) + 
    scale_fill_brewer(palette = "Set2")+
    scale_color_brewer(palette = "Set2")+
    ylim(c(0,100)) + ylab(paste0('% of ', cell_name)) +
    theme_bw() + 
    theme(axis.text.x = element_blank(), legend.position = 'None', axis.title=element_text(size=18),
          axis.text=element_text(size=16, colour="black"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.background = element_blank(), aspect.ratio = 1.6/1.3,
          axis.ticks.length=unit(.25, "cm")) + xlab('') +
    annotate(geom = 'segment', x = 1, xend = 2, y = 70, yend = 70, color = subset(compare_region, variable == cell_type)$line_color[1])+
    geom_text(x = 1.5, y = 73, label = subset(compare_region, variable == cell_type)$p_adj_star[1]) +
    annotate(geom = 'segment', x = 1, xend = 3, y = 85, yend = 85, color = subset(compare_region, variable == cell_type)$line_color[2])+
    geom_text(x = 2, y = 88, label = subset(compare_region, variable == cell_type)$p_adj_star[2]) +
    annotate(geom = 'segment', x = 2, xend = 3, y = 100, yend = 100, color = subset(compare_region, variable == cell_type)$line_color[3])+
    geom_text(x = 2.5, y = 103, label = subset(compare_region, variable == cell_type)$p_adj_star[3])
  
  return(p)
}

plt_sur_by_smokingClass <- function(LUSCmean, class_index, title){
  
  df <- subset(LUSCmean, smokingClass == class_index)
  cut <- surv_cutpoint(df, time = "time", event = "event",
                       variables = c("fi"))
  df_cut <- surv_categorize(cut)
  
  fit <- survfit(Surv(time, event)~fi, data = df_cut)
  p<-ggsurvplot(fit, data = df_cut,
                ylab = "Overall Survival",
                xlab = "Time (Days)",
                break.time.by = 500, xlim = c(0,3000),
                conf.int = F,
                risk.table = T, risk.table.y.text.col = T, risk.table.y.text = F, risk.table.height = 0.3, 
                risk.table.title = "No. at Risk", risk.table.fontsize = 4, 
                pval = T, pval.coord = c(0.1, 0.1,0.1,0.1),
                legend = c(0.85, 0.8), title = title)
  return(p) 
}

cus_theme <-theme(plot.title=element_text(size=16),
                  axis.text=element_text(size=14),
                  axis.title=element_text(size=16),
                  legend.title = element_text(size = 14),
                  legend.text = element_text(size = 14))

#############.
#Data
#############.

load('data/validation_cohort.RData')
load("./data/TCGA_LUAD-LUSC_compath.RData")

lym_list <- c("cd4","cd8","foxp3","cd20","cd20cxcr5","cd79bCoexp")

#Data dictionary ####
#############.
#1# TCGA_LUAD-LUSC_compath.RData: main analysis presented in the paper 
#includes immune hotspot, %lymphocytes, %cancer using tcga luad and lusc
#H&E sections

#LUSCmean: TCGA LUSC cohort with complete survival data (n=462 for 462 patients), including survival data, smoking history, demographic information, disease stages, and hotspot scores.
#stromal_per, lym_per, tumour_per: percentages of stromal cells, lymphocytes and tumour cells in total cells automatically detected in H&E slides.
#fi: immune hotspot score, fraction of cancer-immune hotspot overlaps in immune hotspots.
#fc: cancer hotspot score, fraction of cancer-immune hotspot overlaps in cancer hotspots.
#ft: cancer-immune hotspot score, fraction of cancer-immune hotspot overlaps in the whole tissue.

#LUADmean: TCGA LUAD cohort with complete survival data (n=473 for 473 patients)

#2#. TCGA_LUAD-LUSC_path_EST_ABS_BoLi_DAV.RData: corresponding scores for the same 
#LUAD and LUSC cases with the following datasets:
#TCGA pathology manual scores, 
#ESTIMATE, 
#Davoli, 
#ASBOLUTE, 
#BoLi
#Danaher,
#CIBERSORT. 

#3# TCGA_LUSC-biolinksGeneExp_immuneHotspot.RData: summary of gene expression analysis/enrichment 
#for TCGA LUSC up/down regulated immune genes~immune hotspot groups. 
#alternatively, to reproduce this df, code is provided to use TCGAbiolinks and DESeq2 

#4# TCGA_LUAD-LUSC_compath.RData: main analysis presented in the paper 
#includes immune hotspot, %lymphocytes, %cancer using tcga luad and lusc
#H&E sections 

#5# TCGA_LUAD-LUSC_path_EST_ABS_BoLi_DAV.RData: corresponding scores for the same 
#LUAD and LUSC cases with the following datasets:
#TCGA pathology manual scores, 
#ESTIMATE, 
#Davoli, 
#ASBOLUTE, 
#BoLi
#Danaher,
#CIBERSORT. 

#6# TCGA_LUSC-biolinksGeneExp_immuneHotspot.RData: summary of gene expression analysis/enrichment 
#for TCGA LUSC up/down regulated immune genes~immune hotspot groups. 
#alternatively, to reproduce this df, code is provided to use TCGAbiolinks and DESeq2 

##TRACERx cohort
#7# tx100.RData: computational pathology results from the TRACERx100 diagnostic cohort 

#Validation cohort

#8# - sum_all_no_join: Validation LUSC cohort with serial mIHC sections (n=30 for 10 patients)
#x.s, y.s: coordinates of bottom-right corner of 50x50 μm^2 grids in the slide.
#cell.count.c: number of cancer cells derived from the registered H&E slide.
#cell.count.l: number of lymphocytes derived from the registered H&E slide.
#Hotspots: Cancer, Cancer-immune, Immune hotspots
#Region: TLS: tertiary lymphoid structures, LAG: non-TLS lymphoid aggregates, UD: other regions
#T cell subsets include: cd4: CD4+FOXP3- T cells; cd8: CD8+ T cells; foxp3: CD4+FOXP3+ T cells
#B cell subsets include: cd20: CD20+CXCR5- B cells; cd20cxcr5: CD20+CXCR5+ B cells; cd79bCoexp: CD79b+ B cells; p40: P40+ cancer cells
#uc, hem, cxcr5: stain-negative cells on T and B cell panel and CXCR5 single positive cells 

#9# - sum_all_2: same content as sum_all_no_join except that grids are subgrouped into labelled TLSs/lymphoid aggregates.
#join_idx_r: label of mannual annotated TLSs/lymphoid aggregates (n=130)

#10# - mori_all: SCR (Morisita index of CD8+ and CD4+FOXP3+ T cells) in hotspots.

#11# - df_by_region: number of lymphocytes in 130 mannually labelled TLSs/lymphoid aggregates

#12# - rtre: Target Registration Error between registered and original slides.
#sum_d: sum of distances between benchmarks on the registered and fixed slide.
#mean_rTRE, SD_rTRE: mean and standard diviation of the Target Registration Error representing which represents the distance between original and transformed points normalized by the diagonal length of an image.

#############.
#Fig. 1b ####
#############.
res.cutLUSC <- surv_cutpoint(LUSCmean, time = "time", event = "event", minprop = 0.2,
                             variables = c("fi"))
res.catLUSC <- surv_categorize(res.cutLUSC)
fit <- survfit(Surv(time, event)~fi, data = res.catLUSC)
ggsurvplot(fit, data = res.catLUSC, risk.table = T,  break.time.by = 500, xlim = c(0, 3650),
           pval = T, pval.coord = c(0.1, 0.1))

#############.
#Fig. 1c ####
#############.
res.cutLUSC <- surv_cutpoint(LUSCmean, time = "time", event = "event", minprop = 0.2,
                             variables = c("fi"))
res.catLUSC <- surv_categorize(res.cutLUSC)

LUSCmean$stage_merge = ifelse(LUSCmean$stage %in% c("stage i","stage ia"), "ia", 
                              ifelse(LUSCmean$stage %in% c("stage ib"), "ib",
                                     ifelse(LUSCmean$stage %in% c("stage ii", "stage iia"), "iia",
                                            ifelse(LUSCmean$stage %in% c("stage ib"), "ib",
                                                   ifelse(LUSCmean$stage %in% c("stage iii","stage iiia"), "iiia",
                                                          ifelse(LUSCmean$stage %in% c("stage iiib", "stage iv"), "iv", "iib"))))))

LUSCmean$stage_merge <- as.factor(LUSCmean$stage_merge)

LUSCmean$fiC <- LUSCmean$fi
LUSCmean$fiC[LUSCmean$fiC > res.cutLUSCd[["fi"]]$estimate[[1]]] <- 1
LUSCmean$fiC[LUSCmean$fiC <= res.cutLUSCd[["fi"]]$estimate[[1]]] <- 0

fita <- coxph(Surv(time, event)~fiC+age+stage_merge+pack_years, data = LUSCmean)

fp_a <- ggforest(fita,data = NULL, main = "", fontsize = 0.7, refLabel = "", noDigits = 2)

#############.
#Fig. 1d ####
#############.
ggscatter(LUSCmean, 
          x = "fi", y = "lym_per", 
          add = "reg.line", 
          xlab = "Immune hotspot", 
          ylab = "Lymphocyte%",
          color = "blue4",
          conf.int = TRUE, size = 2,
          add.params = list(color = "grey50", fill = "azure3"), 
          cor.coef = TRUE, cor.method = "spearman")

#############.
#Fig. 1e ####
#############.
load("./data/TCGA_LUAD-LUSC_compath.RData")
res.cutLUSC <- surv_cutpoint(LUSCmean, time = "time", event = "event", minprop = 0.2,
                             variables = c("fi", "lym_per"))
res.catLUSC <- surv_categorize(res.cutLUSC)
res.catLUSC$fi_lym<-NULL
res.catLUSC$fi_lym[res.catLUSC$fi=="low" & res.catLUSC$lym_per=="high"]<-"Low fi, high lym"
res.catLUSC$fi_lym[is.na(res.catLUSC$fi_lym)]<-"All others"
fit <- survfit(Surv(time, event)~fi_lym, data = res.catLUSC)
ggsurvplot(fit, data = res.catLUSC, risk.table = T, break.time.by = 500, xlim = c(0, 3650),
           pval = T, pval.coord = c(0.1, 0.1))

#############.
#Fig. 2a, b ####
#############.

#the below chunk will re-produce already saved DESq2/TCGA-biolinks analysis made using
#the same fi (immune hotspot) cutoff in TCGA LUSC: 
# start with the means
load("./data/TCGA_LUAD-LUSC_compath.RData")
res.cutLUSC <- surv_cutpoint(LUSCmean, time = "time", event = "event", minprop = 0.2,
                             variables = c("fi"))
res.catLUSC <- surv_categorize(res.cutLUSC)
fit <- survfit(Surv(time, event)~fi, data = res.catLUSC)

#this file contains everything saved using the below START gene exp - END gene exp script
load("./data/TCGA_LUSC-biolinksGeneExp_immuneHotspot.RData")

# Enrichment Analysis EA (TCGAVisualize)
# Gene Ontology (GO) and Pathway enrichment barPlot
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP), 
                        GOBPTab = ansEA$ResBP,
                        GOCCTab = ansEA$ResCC,
                        GOMFTab = ansEA$ResMF,
                        PathTab = ansEA$ResPat,
                        nRGTab = Genelist, 
                        nBar = 10,
                        filename="./TCGAvisualize_EAbarplot_fiLUSC_updated.pdf")
TCGAVisualize_volcano(x = dataDEGs$logFC,
                      y = dataDEGs$FDR,
                      filename = "./volcano_fiLUSC_verySigBcell_FCER2.pdf",
                      x.cut = 3,
                      y.cut = 10^-5,
                      names = rownames(dataDEGs),
                      color = c("black","red","darkgreen"),
                      names.size = 2,
                      xlab = " Gene expression fold change (Log2)",
                      legend = "State",
                      highlight = c("CD19", "CD79B", "MS4A1", "CD20", "CXCR5", "FCER2"),
                      title = "Volcano plot (fi-high vs fi-low)",
                      show.names = "both",
                      width = 10)
TCGAVisualize_volcano(x = dataDEGs$logFC,
                      y = dataDEGs$FDR,
                      filename = "./volcano_fiLUSC_verySigBcell.pdf",
                      x.cut = 3,
                      y.cut = 10^-5,
                      names = rownames(dataDEGs),
                      color = c("black","red","darkgreen"),
                      names.size = 2,
                      xlab = " Gene expression fold change (Log2)",
                      legend = "State",
                      highlight = c("CD19", "CD79B", "MS4A1", "CD20", "CXCR5"),
                      title = "Volcano plot (fi-high vs fi-low)",
                      show.names = "both",
                      width = 10)
TCGAVisualize_volcano(x = dataDEGs$logFC,
                      y = dataDEGs$FDR,
                      filename = "./volcano_fiLUSC_verySigBcell3.pdf",
                      x.cut = 1.3,
                      y.cut = 10^-10,
                      names = rownames(dataDEGs),
                      color = c("black","red","darkgreen"),
                      names.size = 2,
                      xlab = " Gene expression fold change (Log2)",
                      legend = "State",
                      highlight = c("CD19", "CD79B", "MS4A1", "CD20", "CXCR5"),
                      title = "Volcano plot (fi-high vs fi-low)",
                      show.names = "both",
                      width = 10)


#if you instead want to use TCGAbiolinks to do an online query for TCGA gene data
#and perform enrichment analysis - here's the code: 

##START gene exp analysis 
# export disc results
LUSCmean$fiC <- LUSCmean$fi
LUSCmean$fiC[LUSCmean$fiC > res.cutLUSC[["fi"]]$estimate[[1]]] <- "high"
LUSCmean$fiC[LUSCmean$fiC <= res.cutLUSC[["fi"]]$estimate[[1]]] <- "low"
listSamples <- as.character(as.factor(LUSCmean$bcr_patient_barcode))
# Query platform Illumina HiSeq with a list of barcode 
query <- GDCquery(project = "TCGA-LUSC", 
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  experimental.strategy = "RNA-Seq",
                  platform = "Illumina HiSeq",
                  file.type = "results",
                  barcode = listSamples, 
                  legacy = TRUE)
# Download a list of barcodes with platform IlluminaHiSeq_RNASeqV2
GDCdownload(query)
# Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
# rsem.genes.results as values
LUSCRnaseqSE <- GDCprepare(query)
LUSCMatrix <- assay(LUSCRnaseqSE,"raw_count") # or LUSCMatrix <- assay(LUSCRnaseqSE,"raw_count")
# For gene expression if you need to see a boxplot correlation and AAIC plot to define outliers you can run
#LUSCRnaseq_CorOutliers <- TCGAanalyze_Preprocessing(LUSCRnaseqSE)
# normalization of genes
dataNorm <- TCGAanalyze_Normalization(tabDF = LUSCMatrix, geneInfo =  geneInfo)
# quantile filter of genes
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)
#discard normal samples
# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt), 
                                   typesample = c("TP"))
# Diff.expr.analysis (DEA)
dataFiltT = dataFilt[,samplesTP]
#check the duplicated ID before the next step!!!
dataFiltT <- dataFiltT[,-c(23)] 
#colnames(dataFiltT) <- substr(colnames(dataFiltT), 1, 12)
# matching for high/low
col <- data.frame(condition = LUSCmean$fiC, 
                  bcr_patient_barcode = LUSCmean$bcr_patient_barcode)
col$bcr_patient_barcode = as.character(as.factor(col$bcr_patient_barcode))
col <- col[col$bcr_patient_barcode %in% substr(colnames(dataFiltT), 1, 12),]
col$id <- col$bcr_patient_barcode
i2 <- match(col$id, substr(colnames(dataFiltT), 1, 12), nomatch=0)
col$id2[i2] <- as.character(colnames(dataFiltT))[i2]
colHigh <- col[ which(col$condition=='high'),]
colLow <- col[ which(col$condition=='low'),]
samplesH <- colHigh$id2
samplesL <- colLow$id2
# Diff.expr.analysis (DEA)
dataDEGs <- TCGAanalyze_DEA(mat1 = dataFiltT[,samplesL],
                            mat2 = dataFiltT[,samplesH],
                            Cond1type = "lowfi",
                            Cond2type = "highfi",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")
# DEGs table with expression values in normal and tumor samples
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,"lowfi","highfi",
                                          dataFiltT[,samplesL],dataFiltT[,samplesH])

# Enrichment Analysis EA
# Gene Ontology (GO) and Pathway enrichment by DEGs list
Genelist <- rownames(dataDEGsFiltLevel)
system.time(ansEA <- TCGAanalyze_EAcomplete(TFname="TCGA LUSC DEA genes, highfi Vs lowfi",Genelist))
#ansEATest <- TCGAanalyze_EAcomplete(TFname="TCGA LUSC DEA genes, highfi Vs lowfi",RegulonList = rownames(dataDEGs)) 

##END gene exp analysis

#############.

#Fig. 2c  ####
#############.
#loading all the following datasets: 
#Danaher, CIBERSORT, MCP, TIMER
load("./data/TCGA_LUAD-LUSC_compath.RData")
load("./data/TCGA_LUAD-LUSC_path_EST_ABS_BoLi_DAV.RData")

res.cutLUSC <- surv_cutpoint(LUSCmean, time = "time", event = "event", minprop = 0.2,
                             variables = c("fi"))
res.catLUSC <- surv_categorize(res.cutLUSC)
LUSCmean$fi.condition[LUSCmean$fi > res.cutLUSC[["fi"]]$estimate[[1]]] <- "High"
LUSCmean$fi.condition[LUSCmean$fi <= res.cutLUSC[["fi"]]$estimate[[1]]] <- "Low"

LUSCfi <- Reduce(function(x, y) merge(x, y, by ="bcr_patient_barcode", all.x=TRUE), list(LUSCmean, danaher))
my_comparisons <- list( c("High", "Low"))
p1 = ggboxplot(LUSCfi, x = "fi.condition", y = "bcell.score.danaher", title = "Danaher et al", 
               xlab = "Immune hotspot", ylab = "B-cell score",
               color = "fi.condition", palette = c("#EE0000", "#00BFFF"), 
               add = "jitter", border = "white")+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox")+
  theme(text = element_text(size=18))+
  theme(legend.position="")+
  theme(plot.title = element_text(hjust = 0.5))

#MCP
load("./data/TCGA_LUAD-LUSC_compath.RData")
load("./data/TCGA_LUAD-LUSC_path_EST_ABS_BoLi_DAV.RData")
res.cutLUSC <- surv_cutpoint(LUSCmean, time = "time", event = "event", minprop = 0.2,
                             variables = c("fi"))
res.catLUSC <- surv_categorize(res.cutLUSC)

fit <- survfit(Surv(time, event)~fi, data = res.catLUSC)
ggsurvplot(fit, data = res.catLUSC, risk.table = T, 
           pval = T, pval.coord = c(0.1, 0.1))
# export disc results
LUSCmean$fi.condition[LUSCmean$fi > res.cutLUSC[["fi"]]$estimate[[1]]] <- "High"
LUSCmean$fi.condition[LUSCmean$fi <= res.cutLUSC[["fi"]]$estimate[[1]]] <- "Low"

LUSCfi <- Reduce(function(x, y) merge(x, y, by ="bcr_patient_barcode", all.x=TRUE), list(LUSCmean, estLUSC2))

#MCP
my_comparisons <- list( c("High", "Low"))
LUSCfi = LUSCfi[!is.na(LUSCfi$MCPcounter_Blineage) ,]
LUSCfi = LUSCfi[! LUSCfi$MCPcounter_Blineage > 5000 ,] #remove the outlier 
p2 = ggboxplot(LUSCfi, x = "fi.condition", y = "MCPcounter_Blineage", title = "MCPCounter",
               xlab = "Immune hotspot", ylab = "B lineage",
               color = "fi.condition", palette = c("#EE0000", "#00BFFF"), 
               add = "jitter", border = "white")+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox")+
  theme(text = element_text(size=18))+
  theme(legend.position="")+
  theme(plot.title = element_text(hjust = 0.5))

#BoLi/TIMER
load("./data/TCGA_LUAD-LUSC_compath.RData")
load("./data/TCGA_LUAD-LUSC_path_EST_ABS_BoLi_DAV.RData")
res.cutLUSC <- surv_cutpoint(LUSCmean, time = "time", event = "event", minprop = 0.2,
                             variables = c("fi"))
res.catLUSC <- surv_categorize(res.cutLUSC)

fit <- survfit(Surv(time, event)~fi, data = res.catLUSC)
ggsurvplot(fit, data = res.catLUSC, risk.table = T, 
           pval = T, pval.coord = c(0.1, 0.1))
# export disc results
LUSCmean$fi.condition[LUSCmean$fi > res.cutLUSC[["fi"]]$estimate[[1]]] <- "High"
LUSCmean$fi.condition[LUSCmean$fi <= res.cutLUSC[["fi"]]$estimate[[1]]] <- "Low"

LUSCfi <- Reduce(function(x, y) merge(x, y, by ="bcr_patient_barcode", all.x=TRUE), list(LUSCmean, BoLi2))
LUSCfi = LUSCfi[!is.na(LUSCfi$TIMER_B_cell) ,]
LUSCfi = LUSCfi[! LUSCfi$TIMER_B_cell > 2 ,] #remove the outlier 
p3 = ggboxplot(LUSCfi, x = "fi.condition", y = "TIMER_B_cell", title = "TIMER",
               xlab = "Immune hotspot", ylab = "B-cell", outlier.shape = NA,
               color = "fi.condition", palette = c("#EE0000", "#00BFFF"), 
               add = "jitter", border = "white")+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox")+
  theme(text = element_text(size=18))+
  theme(legend.position="")+
  theme(plot.title = element_text(hjust = 0.5))

#CIBER
load("./data/TCGA_LUAD-LUSC_compath.RData")
load("./data/TCGA_LUAD-LUSC_path_EST_ABS_BoLi_DAV.RData")

res.cutLUSC <- surv_cutpoint(LUSCmean, time = "time", event = "event", minprop = 0.2,
                             variables = c("fi"))
res.catLUSC <- surv_categorize(res.cutLUSC)
fit <- survfit(Surv(time, event)~fi, data = res.catLUSC)
ggsurvplot(fit, data = res.catLUSC, risk.table = T, 
           pval = T, pval.coord = c(0.1, 0.1))
LUSCmean$fi.condition[LUSCmean$fi > res.cutLUSC[["fi"]]$estimate[[1]]] <- "High"
LUSCmean$fi.condition[LUSCmean$fi <= res.cutLUSC[["fi"]]$estimate[[1]]] <- "Low"

CIBER <- merge(CIBER, LUSCmean, by.y = "bcr_patient_barcode", by.x = "Input.Sample")

p4 = ggboxplot(CIBER, x = "fi.condition", y = "B.cells.memory", title = "CIBERSORT",
               xlab = "Immune hotspot", ylab = "B-cell memory", outlier.shape = NA,
               color = "fi.condition", palette = c("#EE0000", "#00BFFF"), 
               add = "jitter", border = "white")+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox")+
  theme(text = element_text(size=18))+
  theme(legend.position="")+
  theme(plot.title = element_text(hjust = 0.5))

grid.arrange(p1, p4, p2, p3, nrow=1)

#############.
#Fig. 4b ####
#############.

df_cancer_else_join <- sum_all_no_join
levels(df_cancer_else_join$Hotspots) <- c('Else','Cancer-immune','Else','Immune')

df <- df_cancer_else_join %>%
  group_by(slide, Hotspots) %>%
  summarise(cd8 = sum(cd8, na.rm = TRUE), cd4 = sum(cd4, na.rm = TRUE),
            foxp3 = sum(foxp3, na.rm=TRUE), cd20 = sum(cd20, na.rm=TRUE),
            cd20cxcr5 = sum(cd20cxcr5, na.rm = TRUE), cd79bCoexp = sum(cd79bCoexp, na.rm=TRUE))

df[,setdiff(colnames(df), c('slide','Hotspots'))] <- df[,setdiff(colnames(df), c('slide','Hotspots'))] / rowSums(df[,setdiff(colnames(df), c('slide','Hotspots'))]) * 100

df <- melt(df, id = c('slide','Hotspots'))

df$Hotspots <- factor(df$Hotspots , levels = c('Cancer-immune','Immune','Else'))

#statistic
compare_hs <- compare_means(value ~ Hotspots, p.adjust.method = "BH", method='wilcox.test', paired = T,
                            group.by = "variable",
                            data =  df) %>%
              mutate(y_pos = rep(c(60, 75, 90),6),p.adj = format.pval(p.adj, digits = 2),
                     p_adj_star = ifelse(p.adj < 0.001, '***', 
                                         ifelse(p.adj <0.01, '**',
                                                ifelse(p.adj < 0.05, '*', '')))) %>%
              mutate(line_color = ifelse(p_adj_star == "", "white", "black"))

plt_hs_cd8 <- plt_hs_celltype(df, compare_hs, cell_type = 'cd8', cell_name = 'CD8+ T cells')
plt_hs_cd4 <- plt_hs_celltype(df, compare_hs, cell_type = 'cd4', cell_name = 'CD4+FOXP3- T cells')
plt_hs_foxp3 <- plt_hs_celltype(df, compare_hs, cell_type = 'foxp3', cell_name = 'CD4+FOXP3+ T cells')
plt_hs_cd20 <- plt_hs_celltype(df, compare_hs, cell_type = 'cd20', cell_name = 'CD20+CXCR5- B cells')
plt_hs_cd20cxcr5 <- plt_hs_celltype(df, compare_hs, cell_type = 'cd20cxcr5', cell_name = 'CD20+CXCR5+ B cells')
plt_hs_cd79b <- plt_hs_celltype(df, compare_hs, cell_type = 'cd79bCoexp', cell_name = 'CD79b+ B cells')

ggarrange(plt_hs_cd8, plt_hs_cd4, plt_hs_foxp3, plt_hs_cd20, plt_hs_cd20cxcr5, plt_hs_cd79b, ncol = 3, nrow = 2)
#############.

#Fig. 4c ####
#############.

df_by_region_r <- subset(df_by_region, Region !='UD')
df_by_region_r$Region <- as.factor(df_by_region_r$Region)
df_by_region_r$Hotspots <- as.factor(df_by_region_r$Hotspots)
levels(df_by_region_r$Hotspots) <- c('Cancer-immune', 'Cancer-immune', 'Immune','Immune')
df_by_region_r <- rbind(df_by_region_r, subset(df_by_region, Region =='UD' & Hotspots %in% c('Cancer-immune','Immune')))

df_by_region_r$cd8_foxp3 <- df_by_region_r$cd8 / df_by_region_r$foxp3
df_by_region_r[,c(5:10)] <- df_by_region_r[,c(5:10)] / rowSums(df_by_region_r[,c(5:10)]) * 100

df_m <- melt(df_by_region_r, id.vars = c('slide', 'Region', 'X', 'Hotspots'))

df_m$Region <- factor(df_m$Region, levels = c('TLS','LAG','UD'))

compare_region <- compare_means(value ~ Region, p.adjust.method = "BH", method='wilcox.test', paired = F,
                            group.by = "variable",
                            data =  subset(df_m, Hotspots %in% c('Immune'))) %>%
                  mutate(y_pos =rep(c(80, 95, 110),7),p.adj = format.pval(p.adj, digits = 2),
                         p_adj_star = ifelse(p.adj <0.001, '***',
                                                    ifelse(p.adj < 0.01, '**', 
                                                           ifelse(p.adj < 0.05, '*', ifelse(p.adj < 0.0001, '', '****'))))) %>%
                  mutate(line_color = ifelse(p_adj_star == "", "white", "black"))

plt_region_cd8 <- plt_region_celltype(df, compare_region, cell_type = 'cd8', cell_name = 'CD8+ T cells')
plt_region_cd4 <- plt_region_celltype(df, compare_region, cell_type = 'cd4', cell_name = 'CD4+FOXP3- T cells')
plt_region_foxp3 <- plt_region_celltype(df, compare_region, cell_type = 'foxp3', cell_name = 'CD4+FOXP3+ T cells')
plt_region_cd20 <- plt_region_celltype(df, compare_region, cell_type = 'cd20', cell_name = 'CD20+CXCR5- B cells')
plt_region_cd20cxcr5 <- plt_region_celltype(df, compare_region, cell_type = 'cd20cxcr5', cell_name = 'CD20+CXCR5+ B cells')
plt_region_cd79b <- plt_region_celltype(df, compare_region, cell_type = 'cd79bCoexp', cell_name = 'CD79b+ B cells')

ggarrange(plt_region_cd8, plt_region_cd4, plt_region_foxp3, plt_region_cd20, plt_region_cd20cxcr5, plt_region_cd79b, nrow = 2, ncol = 3)
#############.

#Fig. 4d ####
#############.

df <- sum_all_2 %>% group_by(slide, Region, Hotspots, join_idx_r) %>% summarise(n = n()) %>% ungroup() %>% 
      group_by(slide, join_idx_r) %>% slice(which.max(n)) %>% ungroup() %>% 
      group_by(Region, Hotspots) %>% summarise(n_region = n())
#TLS
tls_n <- data.frame(c(df[df[,"Region"] == "TLS", "n_region"]), c(df[df[,"Region"] == "TLS", "Hotspots"]))
tls_n$labs <- paste0(tls_n$n_region, "/", sum(tls_n$n_region))
ggpie(tls_n, "n_region", label = 'labs', lab.pos = "in", lab.font = "white",
      fill = "Hotspots", color = "Gray",
      palette = 'Set2')

#LAG
lag_n <- data.frame(c(df[df[,"Region"] == "LAG", "n_region"]), c(df[df[,"Region"] == "LAG", "Hotspots"]))
lag_n$labs <- paste0(lag_n$n_region, "/", sum(lag_n$n_region))
ggpie(lag_n, "n_region", label = 'labs', lab.pos = "in", lab.font = "white",
      fill = "Hotspots", color = "Gray",
      palette = c("#8DA0CB", "#66C2A5",  "#E78AC3", "#FC8D62"))
#############.

#Fig. 4e ####
#############.

df_region_ratio <- sum_all_no_join %>% group_by(Hotspots, Region, slide) %>% summarise(n = n())

df_region_ratio <- subset(df_region_ratio, Hotspots %in% c('Cancer-immune','Immune'))

df_region_ratio_sum <- df_region_ratio %>% group_by(slide, Hotspots) %>% summarise(sum_n = sum(n))
df_region_ratio <- merge(df_region_ratio, df_region_ratio_sum, by = c('slide', 'Hotspots'))
df_region_ratio$ratio <- df_region_ratio$n / df_region_ratio$sum_n * 100

df_region_ratio <- tidyr::complete(df_region_ratio, slide, Region, Hotspots, fill = list(ratio=0))

#TLS
df_p <- compare_means(ratio ~ Hotspots, p.adjust.method = "BH", method='wilcox.test', paired = T, 
                      data = subset(df_region_ratio, Region == 'TLS' & Hotspots %in% c('Cancer-immune','Immune'))) 

ggboxplot(subset(df_region_ratio, Region == 'TLS' & Hotspots %in% c('Cancer-immune','Immune')), x = 'Hotspots', y = 'ratio', fill = 'Hotspots', palette = c("#EE0000","#00BFFF")) +
  theme_bw() + ylab('% of area') +
  theme( axis.title=element_text(size=18),
         axis.text=element_text(size=16, colour="black"),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank(),
         axis.line = element_line(colour = "black"),
         panel.background = element_blank(), aspect.ratio = 1.6/1.3, 
         axis.ticks.length=unit(.25, "cm"), axis.text.x = element_blank(), axis.title.x = element_blank()) +
  geom_text(x=1.5, y = 15, label=paste0("p = ", df_p$p.adj), color = 'black', size=7)

#LAG
df_p <- compare_means(ratio ~ Hotspots, p.adjust.method = "BH", method='wilcox.test', paired = T, 
                          data = subset(df_region_ratio, Region == 'LAG' & Hotspots %in% c('Cancer-immune','Immune'))) 

ggboxplot(subset(df_region_ratio, Region == 'LAG' & Hotspots %in% c('Cancer-immune','Immune')), x = 'Hotspots', y = 'ratio', fill = 'Hotspots', palette = c("#EE0000","#00BFFF")) +
  theme_bw() + ylab('% of area') +
  theme( axis.title=element_text(size=18),
         axis.text=element_text(size=16, colour="black"),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank(),
         axis.line = element_line(colour = "black"),
         panel.background = element_blank(), aspect.ratio = 1.6/1.3, 
         axis.ticks.length=unit(.25, "cm"), axis.text.x = element_blank(), axis.title.x = element_blank()) +
         geom_text(x=1.5, y = 8.5, label=paste0("p = ", df_p$p.adj), color = 'black', size=7)
  

#############.

#Fig. 5b, c ####
#############.
load("./data/tx100.RData")

sp <- ggplot(diagnostic, aes(x=SCR, y=foxp3_per)) +
  geom_point()
sp + geom_density_2d()

diagnostic$SCR_lym <- NULL
diagnostic$SCR_lym[diagnostic$SCR >= quantile(diagnostic$SCR)[3] &
                     diagnostic$lymphocytes_per >= quantile(diagnostic$lymphocytes_per)[3]] <- "high_high"
diagnostic$SCR_lym[diagnostic$SCR >= quantile(diagnostic$SCR)[3] &
                     diagnostic$lymphocytes_per < quantile(diagnostic$lymphocytes_per)[3]] <- "high_low"
diagnostic$SCR_lym[diagnostic$SCR < quantile(diagnostic$SCR)[3] &
                     diagnostic$lymphocytes_per >= quantile(diagnostic$lymphocytes_per)[3]] <- "low_high"
diagnostic$SCR_lym[diagnostic$SCR < quantile(diagnostic$SCR)[3] &
                     diagnostic$lymphocytes_per < quantile(diagnostic$lymphocytes_per)[3]] <- "low_low"
splots <- list()
fit <- survfit(Surv(DFS_time_days, DFS_censor_variable)~SCR_lym, data = diagnostic)
splots[[1]]<-ggsurvplot(fit, conf.int = FALSE,
                        pval = TRUE, pval.size = 5, pval.coord = c(0.2, 0.1),
                        linetype = "solid",
                        #surv.median.line = "hv",
                        legend = "none", legend.title = "", title = "Tx100",
                        #legend.labs = c("Low", "High"),
                        surv.plot.height = 0.7, palette = get_palette(c("#00AFBB", "#E7B800", "#FC4E07"), 4),
                        risk.table = TRUE,
                        risk.table.col = "black", break.time.by = 200,
                        tables.height = 0.25,
                        tables.theme = theme_cleantable(),
                        tables.y.text = TRUE, risk.table.title = "Number at Risk",
                        tables.x.text = "", xlim = c(0, 1400),
                        xlab = "Days to Death or Recurrence", ylab = "Disease-free Survival")
fit <- survfit(Surv(DFS_time_days, DFS_censor_variable)~SCR_lym, data = diagnostic[diagnostic$Histology =="LUAD",])
splots[[2]]<-ggsurvplot(fit, conf.int = FALSE,
                        pval = TRUE, pval.size = 5, pval.coord = c(0.2, 0.1),
                        linetype = "solid",
                        #surv.median.line = "hv",
                        legend = "none", legend.title = "", title = "LUAD (61)",
                        #legend.labs = c("Low", "High"),
                        surv.plot.height = 0.7, palette = get_palette(c("#00AFBB", "#E7B800", "#FC4E07"), 4),
                        risk.table = TRUE,
                        risk.table.col = "black", break.time.by = 200,
                        tables.height = 0.25,
                        tables.theme = theme_cleantable(),
                        tables.y.text = TRUE, risk.table.title = "Number at Risk",
                        tables.x.text = "", xlim = c(0, 1400),
                        xlab = "Days to Death or Recurrence", ylab = "Disease-free Survival")
fit <- survfit(Surv(DFS_time_days, DFS_censor_variable)~SCR_lym, data = diagnostic[diagnostic$Histology =="LUSC",])
splots[[3]]<-ggsurvplot(fit, conf.int = FALSE,
                        pval = TRUE, pval.size = 5, pval.coord = c(0.2, 0.1),
                        linetype = "solid",
                        #surv.median.line = "hv",
                        legend = "none", legend.title = "", title = "LUSC",
                        #legend.labs = c("Low", "High"),
                        surv.plot.height = 0.7, palette = get_palette(c("#00AFBB", "#E7B800", "#FC4E07"), 4),
                        risk.table = TRUE,
                        risk.table.col = "black", break.time.by = 200,
                        tables.height = 0.25,
                        tables.theme = theme_cleantable(),
                        tables.y.text = TRUE, risk.table.title = "Number at Risk",
                        tables.x.text = "", xlim = c(0, 1400),
                        xlab = "Days to Death or Recurrence", ylab = "Disease-free Survival")
arrange_ggsurvplots(splots, print = TRUE,ncol = 2, nrow = 1)

#############.

#Fig 6a ####
#############.

df <- sum_all_no_join %>% group_by(slide, Hotspots) %>% 
  summarise(cd4 = sum(cd4, na.rm = TRUE), foxp3 = sum(foxp3, na.rm=TRUE))

df$perc <- df$cd4 / df$foxp3

ggpaired(subset(df, Hotspots %in% c('Cancer-immune','Immune')), x= 'Hotspots', y = 'perc', 
         id='slide', line.color = "gray", 
         palette = c("#EE0000","#00BFFF"), color = 'Hotspots',
         xlab = 'Hotspots', ylab = expression("CD8"^"+"~"/"~"CD4"^"+"~"FOXP3"^"+")) + 
  stat_compare_means(paired = TRUE, label.y = 12, label.x=1.4, label = "p.format", method.args = list(alternative = "two.sided")) + 
  rotate_x_text(45) + cus_theme
#############.

#Fig 6b  ####
#############.

ggpaired(mori_all, x= 'Hotspots', y = 'SCR', 
         id='slide', line.color = "gray", 
         palette = c("#EE0000","#00BFFF"), color = 'Hotspots',
         xlab = 'Hotspots', ylab = expression("SCR")) + 
  stat_compare_means(paired = TRUE, label.y = 1, label.x=1.4, label = "p.format", method.args = list(alternative = "two.sided")) + 
  rotate_x_text(45) + cus_theme + theme(axis.text.x = element_blank())
#############.

#Fig 6c ####
#############.

df_den <- subset(sum_all_no_join, Hotspots %in% c('Cancer-immune','Immune'))%>% 
  group_by(slide, Hotspots) %>% 
  summarise(cd8 = sum(cd8, na.rm = TRUE), cd4 = sum(cd4, na.rm = TRUE),
            foxp3 = sum(foxp3, na.rm=TRUE), cd20 = sum(cd20, na.rm=TRUE), 
            cd20cxcr5 = sum(cd20cxcr5, na.rm = TRUE), cd79bCoexp = sum(cd79bCoexp, na.rm=TRUE),
            p40 = sum(p40, na.rm = TRUE), n=n())

df_den$cell_count_all <- rowSums(df_den[,lym_list])
df_den[c('cd20','cd20cxcr5','cd79bCoexp')] <- df_den[c('cd20','cd20cxcr5','cd79bCoexp')]/df_den$cell_count_all * 100

df_den <- merge(df_den, mori_all[c('slide','Hotspots','SCR')], by = c('slide','Hotspots'), all.x = T)
allM <- melt(df_den[c('slide','Hotspots', 'cd20','cd20cxcr5','cd79bCoexp','SCR')], id.vars = c('slide','Hotspots'))
allM <- na.omit(allM)
allM <- subset(allM, Hotspots %in% c('Cancer-immune','Immune'))

allM_mori <- subset(allM, variable=="SCR")
colnames(allM_mori)[which(colnames(allM_mori) == 'value')] <- "cd8_foxp3"
allM <- merge(allM, allM_mori[,c('slide','Hotspots','cd8_foxp3')], by = c('slide','Hotspots'))

levels(allM$variable) <- c("CD20+CXCR5- B cells", "CD20+CXCR5+ B cells", "CD79b+ B cells", "Morisita")

ggscatter(subset(allM, variable != 'Morisita'), x = 'value', y = 'cd8_foxp3', add = "reg.line",  
          facet.by = c('Hotspots','variable'), scales = 'free_x',# Add regression line
          color = 'Hotspots',conf.int = TRUE, palette = c("#EE0000","#00BFFF"),                                 # Add confidence interval
          add.params = list(color = "gray40",
                            fill = "lightgray"),
          xlab = 'Cell percentage(%)', ylab= "SCR")+
  stat_cor(method = "pearson", label.y = 1) + ylim(c(0.6,1.05))+
  theme(plot.title=element_text(size=14),
        axis.text=element_text(size=7),
        axis.title=element_text(size=16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))

#############.

#Fig. S1a, b, c + Fig. S3d: corr maps  ####
#############.
load("./data/TCGA_LUAD-LUSC_compath.RData")
#loading all the following datasets: 
#LUAD.TCGAPathScores, estLUAD2, Davoli, absolute2, BoLi2
load("./data/TCGA_LUAD-LUSC_path_EST_ABS_BoLi_DAV.RData")


#merges
LUADtumor <- Reduce(function(x, y) merge(x, y, by ="bcr_patient_barcode", all.x=TRUE), list(LUADmean, LUAD.TCGAPathScores, estLUAD2, absolute2))
LUSCtumor <- Reduce(function(x, y) merge(x, y, by ="bcr_patient_barcode", all.x=TRUE), list(LUSCmean, LUSC.TCGAPathScores, estLUSC2, absolute2))
LUADimmune <- Reduce(function(x, y) merge(x, y, by ="bcr_patient_barcode", all.x=TRUE), list(LUADmean, estLUAD2, BoLi2, Davoli))
LUSCimmune <- Reduce(function(x, y) merge(x, y, by ="bcr_patient_barcode", all.x=TRUE), list(LUSCmean, estLUSC2, BoLi2, Davoli))
LUADstroma <- Reduce(function(x, y) merge(x, y, by ="bcr_patient_barcode", all.x=TRUE), list(LUADmean, LUAD.TCGAPathScores, estLUAD2))
LUSCstroma <- Reduce(function(x, y) merge(x, y, by ="bcr_patient_barcode", all.x=TRUE), list(LUSCmean, LUSC.TCGAPathScores, estLUSC2))
#prepare
rownames(LUADtumor) <- LUADtumor$bcr_patient_barcode
rownames(LUSCtumor) <- LUSCtumor$bcr_patient_barcode
rownames(LUADimmune) <- LUADimmune$bcr_patient_barcode
rownames(LUSCimmune) <- LUSCimmune$bcr_patient_barcode
rownames(LUADstroma) <- LUADstroma$bcr_patient_barcode
rownames(LUSCstroma) <- LUSCstroma$bcr_patient_barcode

LUADtumor <- subset(LUADtumor, select = c(tumour_per, TCGA_pathTumorPer, ARAN_ESTIMATE, ARAN_ABSOLUTE) )
LUSCtumor <- subset(LUSCtumor, select = c(tumour_per, TCGA_pathTumorPer, ARAN_ESTIMATE, ARAN_ABSOLUTE) )

LUSCimmunefiMCPONLY <-subset(LUSCimmune, select = c(fi, MCPcounter_Tcells, MCPcounter_CD8Tcells, MCPcounter_Cytotoxiclymphocytes, MCPcounter_NKcells, MCPcounter_Blineage, 
                                                    MCPcounter_Monocyticlineage, MCPcounter_Myeloiddendriticcells, MCPcounter_Neutrophils, MCPcounter_Endothelialcells))

LUSCimmunefi <-subset(LUSCimmune, select = c(fi, DAVOLI_TotNofMutations.in.exons, DAVOLI_TP53Mutation, DAVOLI_CellCycle.Signature.Score, DAVOLI_SCNA.Level, DAVOLI_Chrom.SCNA.Level, 
                                             DAVOLI_Arm.SCNA.Level, DAVOLI_Focal.SCNA.Level, DAVOLI_Chrom.Arm.SCNA.Level, DAVOLI_SCNA.Level.normalized.by.size, 
                                             MCPcounter_Tcells, MCPcounter_CD8Tcells, MCPcounter_Cytotoxiclymphocytes, MCPcounter_NKcells, MCPcounter_Blineage, 
                                             MCPcounter_Monocyticlineage, MCPcounter_Myeloiddendriticcells, MCPcounter_Neutrophils, MCPcounter_Endothelialcells))

LUADimmune <-subset(LUADimmune, select = c(fi, lym_per, ESTIMATE_ImmuneScore, MCPcounter_immune_sum, TIMER_B_cell, TIMER_T_cell.CD4, TIMER_T_cell.CD8, 
                                           TIMER_Neutrophil, TIMER_Macrophage, TIMER_DC, DAVOLI_TotNofMutations.in.exons, DAVOLI_TP53Mutation, 
                                           DAVOLI_Immune.Signature.Score, DAVOLI_CellCycle.Signature.Score, DAVOLI_SCNA.Level, DAVOLI_Chrom.SCNA.Level, 
                                           DAVOLI_Arm.SCNA.Level, DAVOLI_Focal.SCNA.Level, DAVOLI_Chrom.Arm.SCNA.Level, DAVOLI_SCNA.Level.normalized.by.size))
LUSCimmune <-subset(LUSCimmune, select = c(fi, lym_per, ESTIMATE_ImmuneScore, MCPcounter_immune_sum, TIMER_B_cell, TIMER_T_cell.CD4, TIMER_T_cell.CD8, 
                                           TIMER_Neutrophil, TIMER_Macrophage, TIMER_DC, DAVOLI_TotNofMutations.in.exons, DAVOLI_TP53Mutation, 
                                           DAVOLI_Immune.Signature.Score, DAVOLI_CellCycle.Signature.Score, DAVOLI_SCNA.Level, DAVOLI_Chrom.SCNA.Level, 
                                           DAVOLI_Arm.SCNA.Level, DAVOLI_Focal.SCNA.Level, DAVOLI_Chrom.Arm.SCNA.Level, DAVOLI_SCNA.Level.normalized.by.size))


LUADstroma <-subset(LUADstroma, select = c(stromal_per, TCGA_pathStromalPer, ESTIMATE_StromalScore, MCPcounter_stromal_sum))
LUSCstroma <-subset(LUSCstroma, select = c(stromal_per, TCGA_pathStromalPer, ESTIMATE_StromalScore, MCPcounter_stromal_sum))

tumour <- rbind(LUADtumor, LUSCtumor)
stroma <- rbind(LUADstroma, LUSCstroma)
immune <- rbind(LUADimmune, LUSCimmune)
immune <- subset(immune, select = -c(fi, DAVOLI_TotNofMutations.in.exons, DAVOLI_TP53Mutation, DAVOLI_SCNA.Level, DAVOLI_SCNA.Level, DAVOLI_Chrom.SCNA.Level, 
                                     DAVOLI_Arm.SCNA.Level, DAVOLI_Focal.SCNA.Level, DAVOLI_Chrom.Arm.SCNA.Level, DAVOLI_SCNA.Level.normalized.by.size, DAVOLI_CellCycle.Signature.Score))

tumour <- na.omit(tumour)
stroma <- na.omit(stroma)
immune <- na.omit(immune)
LUADtumor <- na.omit(LUADtumor)
LUSCtumor <- na.omit(LUSCtumor)
LUADimmune <- na.omit(LUADimmune)
LUSCimmune <- na.omit(LUSCimmune)
LUADstroma <- na.omit(LUADstroma)
LUSCstroma <- na.omit(LUSCstroma)
LUSCimmunefi <- na.omit(LUSCimmunefi)
LUSCimmunefiMCPONLY <- na.omit(LUSCimmunefiMCPONLY)

### Corr maps
corTumor <- round(cor(tumour),2)
corStroma <- round(cor(stroma),2)
corImmune <- round(cor(immune),2)
corLUADtumor <- round(cor(LUADtumor),2)
corLUSCtumor <- round(cor(LUSCtumor),2)
corLUADimmune <- round(cor(LUADimmune),2)
corLUSCimmune <- round(cor(LUSCimmune),2)
corLUADstroma <- round(cor(LUADstroma),2)
corLUSCstroma <- round(cor(LUSCstroma),2)
corLUSCimmunefi <- round(cor(LUSCimmunefi),2)
corLUSCimmunefiMCPONLY <- round(cor(LUSCimmunefiMCPONLY),2)

ggcorrplot(corTumor, method = "circle", title = "TCGA NSCLC - cancer cell, n=401")
ggcorrplot(corStroma, method = "circle", title = "TCGA NSCLC - stroma cell, n=880")
ggcorrplot(corImmune, method = "circle", title = "TCGA NSCLC - immune cell, n=689")
ggcorrplot(corLUSCimmunefiMCPONLY, method = "circle", title = "TCGA LUSC - immune hotspot, n=463")

#############.

##Fig. S3a, b, c ####
#############.
#same as above we need the same data, or tcga query to download again - as shown above in Fig 2
load("./data/TCGA_LUAD-LUSC_compath.RData")
load("./data/TCGA_LUSC-biolinksGeneExp_immuneHotspot.RData")

#exhaustedB signature 
grep("CD27", rownames(dataFiltT), value = TRUE)
exhaustedB <- as.data.frame( dataFiltT[c("CD19", "CD4", "CD69", "CD27", "MS4A1", "FOXP3", 
                                         "CD8A", "CD8B" ,"CD80" ,"CD86", "CD81" ,"CD82", "CD84", "CD83", "CR2" ),])
exhaustedB <- as.data.frame(t(exhaustedB))
exhaustedB$id2 <- rownames(exhaustedB)

exhaustedB <- merge(exhaustedB, col, by = "id2")
#take out the outlier in row 381 TCGA-85-A4PA-01A-11R-A24Z-07
#exhaustedB <- exhaustedB[-c(381),]

ggscatter(exhaustedB, x = "FOXP3", y = "CD27", color = "condition",
          add = "reg.line", conf.int = TRUE, size = 2,
          add.params = list(color = "blue",
                            fill = "lightgray"), title = "TCGA LUSC (n = 466)",
          cor.coef = TRUE, cor.method = "spearman")
my_comparisons <- list( c("high", "low"))
ggboxplot(exhaustedB, x = "condition", y = "FOXP3", add = "jitter", color = "condition", title = "TCGA LUSC (n=466)",
          palette = c("#EE0000", "#00BFFF")) +stat_compare_means(comparisons = my_comparisons)
ggboxplot(exhaustedB, x = "condition", y = "CD27", add = "jitter", color = "condition", title = "TCGA LUSC (n=466)",
          palette = c("#EE0000", "#00BFFF")) +stat_compare_means(comparisons = my_comparisons)

#take out the outlier in row 381 TCGA-85-A4PA-01A-11R-A24Z-07
exhaustedBOutlier <- exhaustedB[-c(381),]
exhaustedBOutlier <- exhaustedBOutlier[order(exhaustedBOutlier$condition),]  #order for the heatmap!
annot <- subset(exhaustedBOutlier, select = c(id2, condition))
rownames(annot) <- annot$id2
annot <- subset(annot, select = -c(id2))

mapBcell <- subset(exhaustedBOutlier, select = -c(bcr_patient_barcode, id, condition, FOXP3CD8A.Pheno, meanExhaustedB.cell, 
                                                  CD8B, CD80, CD86, CD81, CD82, CD84, CD83))
corBcell <- subset(mapBcell, select = -c(id2))
rownames(mapBcell) <- mapBcell$id2
mapBcell <- subset(mapBcell, select = -c(id2))
mapBcell <- t(mapBcell)
condition = c("#EE0000", "#00BFFF")
names(condition) = c("high", "low")
ann_colors = list(condition = condition)

mapBcelllog <- log10(mapBcell+1)
pheatmap(mapBcelllog, annotation = annot,annotation_legend = T, annotation_colors = ann_colors, show_colnames= F, 
         gaps_col =  3, cluster_rows=T, cluster_cols=T, drop_levels = F, border_color = NA)


grid.ls(grid.force()) # "col_annotation" looks like it's the one to edit
grid.gedit("col_annotation", gp = gpar(col="grey70"))
matcorBcell <- round(cor(corBcell),2)

ggcorrplot(matcorBcell, method = "square", hc.order = TRUE, type = "upper", lab=TRUE,
           outline.col = "white", title = "TCGA LUSC, n=462")

#############.

#Fig S4a ####
#############.

rtre$mean_rTRE <- rtre$mean_rTRE*100
ggdensity(rtre, x = 'mean_rTRE', rug = TRUE, add = 'mean', xlab = 'Mean Target Registration Error (%)', ylab = 'Density', color = "tomato2") + cus_theme

#############.
#Fig S4b ####
#############.

sum_df <- data.frame(slide  = sum_all_no_join$slide, hs = sum_all_no_join$Hotspots, p40=sum_all_no_join$p40, 
                     cell.count.c = sum_all_no_join$cell.count.c, cell.count.l =sum_all_no_join$cell.count.l,
                     lym_ihc = rowSums(sum_all_no_join[,lym_list]))

df <- subset(sum_df, hs %in% c('Cancer-immune','Immune'))  %>% 
  group_by(slide,hs) %>% 
  summarise(p40 = mean(p40,na.rm=TRUE), cell.count.c = mean(cell.count.c,na.rm=TRUE), 
            cell.count.l = mean(cell.count.l, na.rm=TRUE),
            lym_ihc = mean(lym_ihc, na.rm=TRUE))

df[,setdiff(colnames(df), c("slide","hs"))] <- df[,setdiff(colnames(df), c("slide","hs"))] / (50*50) * 10^6
#levels(df$hs) <- c('Cancer','Cancer-immune','Else','Immune')
colnames(df)[which(colnames(df)=='hs')] <- 'Hotspots'

#cancer
scatter_cancer <- ggscatter(df, x = 'cell.count.c', y = 'p40', add = "reg.line",                
          conf.int = TRUE, color='Hotspots',  palette = c("#EE0000","#00BFFF"),              
          add.params = list(color = "gray40",
                            fill = "lightgray"), 
          ylab = expression("P40"^"+"~"cell / "~mm^2),
          xlab = expression("H&E-based cancer cell / "~mm^2)) + stat_cor(method = "pearson") + cus_theme + coord_fixed(ratio = 1)

#lymphocyte
scatter_lym <- ggscatter(df, x = 'cell.count.l', y = 'lym_ihc', add = "reg.line",                     
          conf.int = TRUE,  color = 'Hotspots', palette = c("#EE0000","#00BFFF"),              
          add.params = list(color = "gray40",
                            fill = "lightgray"), 
          xlab = expression("H&E-based lymphocyte cell / "~mm^2),
          ylab = expression("IHC-based lymphocyte cell / "~mm^2)) + stat_cor(method = "pearson") + cus_theme + coord_fixed(ratio = 0.6)

ggarrange(scatter_cancer, scatter_lym, ncol = 2)
#############.

#Fig S4c ####
#############.
#density of cell types in the whole slide

df <- aggregate(. ~ x.s+y.s+slide, sum_all_no_join, sum)
df <- df %>% group_by(slide) %>% summarise(cd8 = sum(cd8, na.rm = TRUE), cd4 = sum(cd4, na.rm = TRUE),
                                           foxp3 = sum(foxp3, na.rm=TRUE), cd20 = sum(cd20, na.rm=TRUE), 
                                           cd20cxcr5 = sum(cd20cxcr5, na.rm = TRUE), cd79bCoexp = sum(cd79bCoexp, na.rm=TRUE),
                                           p40 = sum(p40, na.rm = TRUE), n=n())
df[,lym_list] <- df[,lym_list]/ df$n /(50*50) * 10^6

df_den <- melt(df[,c('slide', lym_list)], id.vars = c('slide'))

df_den <- tidyr::complete(df_den, slide, variable, fill = list(value=NA))
df_den$variable <- factor(df_den$variable, levels = c('cd8','cd4','foxp3','cd20','cd20cxcr5','cd79bCoexp'))

df_p_den <- compare_means(value ~ variable, p.adjust.method = "BH", method='wilcox.test', paired = T, 
                          #group.by =  c("slide"),
                          data =  df_den) #%>% mutate(y_pos =rep(c(6000, 7000, 8000),6))
df_p_compare <- list(c('cd8','cd4'), c('cd4','foxp3'), c('foxp3','cd20'), 
                     c('cd20','cd20cxcr5'),c('cd20cxcr5','cd79bCoexp'), c('cd20','cd79bCoexp'),
                     c('cd8','foxp3'))

df_p_den <- compare_means(value ~ variable, p.adjust.method = "BH", method='wilcox.test', paired = T, 
                          data = df_den) 

select_row <- c()
for (pair in df_p_compare){
  select_row <- c(select_row, which((df_p_den$group1 == pair[1] & df_p_den$group2 == pair[2]) | (df_p_den$group1 == pair[2] && df_p_den$group2 == pair[1])))
}

df_p_den <- df_p_den[select_row,] %>%
              mutate(y_pos = c(1500, 1500, 1200, 800, 800, 1000, 1700), p.adj = format.pval(p.adj, digits = 1),
                     p_adj_star = ifelse(p.adj < 0.001, '***', 
                                         ifelse(p.adj <0.01, '**',
                                                ifelse(p.adj < 0.05, '*', 'ns'))))
ggboxplot(df_den, 
          x= 'variable', y = 'value', 
          id='slide',
          add = 'jitter', add.params = list(size=0.8), ylim = c(0,2000),
          xlab = 'Immune cell subsets', ylab=expression("Cell /"~mm^2)) +
  #stat_compare_means(paired = TRUE, comparisons = list(c("LAG","TLS"), c('LAG','UD'), c('TLS','UD')), label = "p.format") +
  # stat_compare_means(paired = TRUE,
  #                    label = "p.signif", comparisons = df_p_compare,  label.y = c(1500, 1500, 1200, 800, 800, 1000, 1700)) +
  geom_signif(data=df_p_den,
              aes(xmin = group1, xmax = group2, annotations = p_adj_star, y_position = y_pos),
              manual= TRUE) +
  scale_x_discrete(labels = c(expression("CD8"^"+"),
                   expression("CD4"^"+"~"FOXP3"^"-"),
                   expression("CD4"^"+"~"FOXP3"^"+"),
                   expression("CD20"^"+"~"CXCR5"^"-"),
                   expression("CD20"^"+"~"CXCR5"^"+"),
                   expression("CD79b"^"+")))+
  rotate_x_text(45) + cus_theme
#############.

#Fig S5 ####
#############.

#LUAD
p_sur_luad_1 <- plt_sur_by_smokingClass(LUADmean, "1", "Current reformed smoker for <= 15 years")
p_sur_luad_2 <- plt_sur_by_smokingClass(LUADmean, "2", "Current reformed smoker for > 15 years")
p_sur_luad_3 <- plt_sur_by_smokingClass(LUADmean, "3", "Current smoker")
p_sur_luad_4 <- plt_sur_by_smokingClass(LUADmean, "4", "Lifelong non-smoker")

LUADmean_1 <- subset(LUADmean, smokingClass %in% c('1','2','3','4'))

LUADmean_1$stage_merge = ifelse(LUADmean_1$stage %in% c("stage i","stage ia"), "ia", 
                              ifelse(LUADmean_1$stage %in% c("stage ib"), "ib",
                                     ifelse(LUADmean_1$stage %in% c("stage ii", "stage iia"), "iia",
                                            ifelse(LUADmean_1$stage %in% c("stage ib"), "ib",
                                                   ifelse(LUADmean_1$stage %in% c("stage iii","stage iiia"), "iiia",
                                                          ifelse(LUADmean_1$stage %in% c("stage iiib", "stage iv"), "iv", "iib"))))))

res.cutLUADd <- surv_cutpoint(LUADmean_1, time = "time", event = "event",
                              variables = c("fi"))
res.catLUADd <- surv_categorize(res.cutLUADd)

LUADmean_1$fiC <- LUADmean_1$fi
LUADmean_1$fiC[LUADmean_1$fiC > res.cutLUADd[["fi"]]$estimate[[1]]] <- 1
LUADmean_1$fiC[LUADmean_1$fiC <= res.cutLUADd[["fi"]]$estimate[[1]]] <- 0

fita <- coxph(Surv(time, event)~fiC+stage_merge+smokingClass, data = LUADmean_1)
ggforest(fita, data = NULL, main = "", fontsize = 0.7, refLabel = "", noDigits = 2)

#LUSC
p_sur_lusc_1 <- plt_sur_by_smokingClass(LUSCmean, "1", "Current reformed smoker for <= 15 years")
p_sur_lusc_2 <- plt_sur_by_smokingClass(LUSCmean, "2", "Current reformed smoker for > 15 years")
p_sur_lusc_3 <- plt_sur_by_smokingClass(LUSCmean, "3", "Current smoker")
p_sur_lusc_4 <- plt_sur_by_smokingClass(LUSCmean, "4", "Lifelong non-smoker")

LUSCmean_1 <- subset(LUSCmean, smokingClass %in% c('1','2','3','4'))

LUSCmean_1$stage_merge = ifelse(LUSCmean_1$stage %in% c("stage i","stage ia"), "ia", 
                                ifelse(LUSCmean_1$stage %in% c("stage ib"), "ib",
                                       ifelse(LUSCmean_1$stage %in% c("stage ii", "stage iia"), "iia",
                                              ifelse(LUSCmean_1$stage %in% c("stage ib"), "ib",
                                                     ifelse(LUSCmean_1$stage %in% c("stage iii","stage iiia"), "iiia",
                                                            ifelse(LUSCmean_1$stage %in% c("stage iiib", "stage iv"), "iv", "iib"))))))

res.cutLUSCd <- surv_cutpoint(LUSCmean_1, time = "time", event = "event",
                              variables = c("fi"))
res.catLUSCd <- surv_categorize(res.cutLUSCd)

LUSCmean_1$fiC <- LUSCmean_1$fi
LUSCmean_1$fiC[LUSCmean_1$fiC > res.cutLUSCd[["fi"]]$estimate[[1]]] <- 1
LUSCmean_1$fiC[LUSCmean_1$fiC <= res.cutLUSCd[["fi"]]$estimate[[1]]] <- 0

fita <- coxph(Surv(time, event)~fiC+stage_merge+smokingClass, data = LUSCmean_1)
ggforest(fita, data = NULL, main = "", fontsize = 0.7, refLabel = "", noDigits = 2)

#############.
