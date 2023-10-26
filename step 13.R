## HPA analysis code

rm(list=ls())
setwd('D:/Rbook')
#BiocManager::install('BiocStyle')
library(BiocStyle)
#BiocManager::install('HPAanalyze')
library(HPAanalyze)
library(dplyr)
library(tibble)
library(readr)
library(tidyr)
dir.create("./step1_img/")
gene='bmp4'
tissue='Breast'

getOption('timeout')
options(timeout=10000)
hpa_target_gene<-hpaXmlGet(gene)
#
hpa_target_gene_fig_url<-hpaXmlTissueExpr(hpa_target_gene)
hpa_target_gene_fig_url_1<-as.data.frame(hpa_target_gene_fig_url[[1]])
hpa_target_gene_fig_url_1[1:6,1:18]

hpa_target_gene_fig_url_2<-as.data.frame(hpa_target_gene_fig_url[[2]])
hpa_target_gene_fig_url_2[1:6,1:18]
#
hpa_target_gene_fig_url_tissue<-hpa_target_gene_fig_url_1[hpa_target_gene_fig_url_1$tissueDescription2==tissue,]
hpa_target_gene_fig_url_tissue<-hpa_target_gene_fig_url_2[hpa_target_gene_fig_url_2$tissueDescription2==tissue,]
#
picDir <- paste('./step1_img/',gene, tissue,"IHC-2/", sep = "_")
if (!dir.exists(picDir)) {
  dir.create(picDir)
}


for (i in 1:nrow(hpa_target_gene_fig_url_tissue)) {
  file_url<-hpa_target_gene_fig_url_tissue$imageUrl[i]
  file_dir<-paste(picDir,gene,tissue,hpa_target_gene_fig_url_tissue$patientId[i],hpa_target_gene_fig_url_tissue$tissueDescription1[i],hpa_target_gene_fig_url_tissue$tissueDescription2[i],".tiff",sep = "_")
  download.file(url = file_url,destfile = file_dir,mode = "wb")
}
write.csv(hpa_target_gene_fig_url_tissue,paste(picDir,gene,"IHC-2_result_tab.csv",sep = "_"))