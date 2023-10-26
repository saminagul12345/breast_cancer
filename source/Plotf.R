

suppressMessages(library(survival))
suppressMessages(library(survminer))
suppressMessages(library(ggplotify))
suppressMessages(library(cowplot))
suppressMessages(library(Hmisc))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(gridExtra))


factorf <- function(fac_mix,factGene,intype){
  
  pdf(file = paste0('fig/step_7_',intype,'_fator_p3.pdf'),width = 6,height = heig)
  print(p3)
  dev.off()
  
}