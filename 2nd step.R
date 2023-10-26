 # 2nd step
rm(list = ls())
options(stringsAsFactors = F)

suppressMessages(library(pacman))
source('source/LoadT_Datf.R')
source('source/oclr.R')
source('source/immunef.R')
source('source/fpkm2Tpm.R')


## Load data: 
TCGA_EXP <- Readdat(inpath = 'input/TCGA-BRCA.htseq_counts.tsv.gz',
                    idpath = 'input/gencode.v22.annotation.gene.probeMap')
TCGA_clinc <- ReadClinc(inpath = 'input/tnbc_exp_cli.Rdata')

#save separately OCLR and ESTIMATE read in file:
TCGA_oclr <- cbind(gene_id = rownames(TCGA_EXP),TCGA_EXP)
write.table(TCGA_oclr,file = 'tmp/step_1_TCGA_oclr.txt',sep = '\t',quote = F)
write.table(TCGA_EXP,file = 'tmp/step_1_TCGA_EXP.txt',sep = '\t',quote = F)


## Clinical data: Extract useful clinical information
TCGA_clinc$OS.time <- as.numeric(signif(TCGA_clinc$days_to_last_follow_up.diagnoses/30,3))
TCGA_clinc$OS <- ifelse(TCGA_clinc$vital_status.demographic == 'Alive',0,
                        ifelse(TCGA_clinc$vital_status.demographic == 'Dead',1,'NA'))
TCGA_clinc$OS <- as.numeric(TCGA_clinc$OS)
TCGA_clinc$Status <- ifelse(TCGA_clinc$OS == 0,'Alive',ifelse(TCGA_clinc$OS == 1,'Dead','NA'))
TCGA_clinc$Age <- TCGA_clinc$age_at_initial_pathologic_diagnosis
TCGA_clinc$Age2 <- ifelse(TCGA_clinc$Age > 60,'>60','<=60')
TCGA_clinc$Stage <- str_extract(TCGA_clinc$clinical_stage,pattern = "[IV]+")
TCGA_clinc$Stage2 <- ifelse(TCGA_clinc$Stage %in% c('I','II'),'I/II','III/IV')



## Calculate the fraction:StenScore
mRNAsi <- oclrpd(sigpath = 'input/pcbc-stemsig.tsv',
                 datpath = 'tmp/step_1_TCGA_oclr.txt')

## Calculate the immune score: ESTIMATE
ESTIscore <- Immestif(datpath = 'tmp/step_1_TCGA_EXP.txt',inplat = 'illumina',insep = 1)

## Organize data
TCGA_clinc <- Reduce(function(x,y) inner_join(x,y,by="submitter",all.x=TRUE),
                     list(TCGA_clinc,mRNAsi,ESTIscore),accumulate =FALSE)

## Save the data
save(TCGA_EXP, file = 'rdata/step_1_TCGA_EXP.Rdata')
save(TCGA_clinc,file = 'rdata/step_1_TCGA_clinc.Rdata')