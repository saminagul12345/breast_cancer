rm(list = ls())
options(stringsAsFactors = F)

suppressMessages(library(dplyr))
source('source/subclusterf.R')

## load the data
load('rdata/step_1_TCGA_EXP.Rdata')
load('rdata/step_4_TCGA_clinc.Rdata')
load('rdata/step_4_DEG_INFO.Rdata')

##  Organize data
DEG_dat <- TCGA_EXP[DEG_INFO[DEG_INFO$State != 'not',]$SYMBOL,]

##  Unsupervised clustering
newcluster(indat = DEG_dat,ink = 10)

##  Add to clinical information
Kcluster <- read.csv('cluster/cluster.k=2.consensusClass.csv',header = F)
names(Kcluster) <- c('submitter','cluster')

TCGA_clinc <- left_join(TCGA_clinc,Kcluster,by = 'submitter')
TCGA_clinc$StemType <- ifelse(TCGA_clinc$cluster == 1,'Steamness Subtype I','Steamness Subtype II')

##  save data
save(TCGA_clinc,file = 'rdata/step_5_TCGA_clinc.Rdata')