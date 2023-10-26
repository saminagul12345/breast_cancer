##first_step_of_project extract TNBC

getwd()

##TCGA data analysis
install.packages("DT")
library(DT)
options(stringsAsFactors = F)
dataPath="/TCGA-BRCA.GDC_phenotype.csv"
dataPath="/home/datahup/mina/data_analysis/R_data/my_Rscript/TNBC2022/tcga/"

setwd(dataPath)

required_pkgs<-c("tidyverse","DESeq2","limma","ggplot2","RColorBrewer",
                 "gplots","amap","BiocParallel")
install.packages("tidyverse")
library(data.table)
install.packages("R.utils")
raw_data <- fread("/opt/disease/TCGA/BRCA/TCGA-BRCA.htseq_counts.tsv.gz", sep = "\t", header = T)
raw_data <- as.data.frame(raw_data)
rownames(raw_data) <- raw_data[,1]
raw_data <- raw_data[,-1]
raw_data <- 2^raw_data-1
raw_data=ceiling(raw_data)
#raw_data[1:4,1:4]
Ensembl_ID = substr(rownames(raw_data), start = 1, stop = 15)
#Ensembl_ID[1:4]
rownames(raw_data)<-Ensembl_ID
#raw_data[1:4,1:4]
save(raw_data, file = 'TCGA-BRCA.htseq_counts.Rdata')


length(raw_data[1,])
raw_data[1:4,1:4]
dim(raw_data)
read


phe3 <- read.delim(file = '/home/datahup/mina/tnbc/TCGA.BRCA.sampleMap_BRCA_clinicalMatrix')
table(phe3$breast_carcinoma_estrogen_receptor_status)
table(phe3$breast_carcinoma_progesterone_receptor_status)
table(phe3$lab_proc_her2_neu_immunohistochemistry_receptor_status)

phe4 <- data.frame(Sample = phe3$sampleID, 
                   ER = phe3$breast_carcinoma_estrogen_receptor_status,
                   PR = phe3$breast_carcinoma_progesterone_receptor_status,
                   HER2 = phe3$lab_proc_her2_neu_immunohistochemistry_receptor_status,
                   last_follwup = phe3$days_to_last_followup,
                   day_dead = phe3$days_to_death,
                   MState = phe3$pathologic_M,
                   NState = phe3$pathologic_N,
                   TState = phe3$pathologic_T,
                   State = phe3$pathologic_stage,
                   OS = phe3$vital_status,
                   age = phe3$age_at_initial_pathologic_diagnosis)
str(phe4)

# phe4$subtype <- apply(phe4, 1, function(x) sum(x =='Negative'))
# table(phe4$subtype)
phe4$subtype <- ifelse(phe4$ER == 'Negative' & phe4$PR == 'Negative' & phe4$HER2 == 'Negative', 'TNBC',
                       ifelse(phe4$ER == 'Positive'  & phe4$HER2 == 'Negative', 'LumA',
                              ifelse(phe4$ER == 'Positive'  & phe4$HER2 == 'Positive', 'LumB',
                                     ifelse(phe4$ER == 'Negative' & phe4$HER2 == 'Positive', 'HER2', 'unknow'))))

table(phe4$subtype)
rownames(phe4) <- phe4$Sample
names(phe4)[1] <- 'Id'
phe0 = subset(phe4,phe4$subtype == 'TNBC')
save(phe0,file = 'tnbc_clinic.Rdata')


load(file = 'tnbc_clinic.Rdata')
exp = raw_data
exp[1:4,1:4]
rownames(phe0)[1:4]
library(stringi)
colnames(exp) = str_sub(colnames(exp), start = 1,end = 15)
exp[1:4,1:4]
exp_tnbc = exp[,colnames(exp) %in% rownames(phe0)]

save(exp_tnbc,phe0,file = 'tnbc_exp_cli.Rdata')


rm(list = ls())
load(file = 'tnbc_exp_cli.Rdata')
load(file ='TCGA-BRCA.htseq_counts.Rdata' )