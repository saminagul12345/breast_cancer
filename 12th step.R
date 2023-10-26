## 12th step scRNA

rm(list = ls())
setwd('/home/datahup/pjy/BRCA/tang/mina/')


library(Seurat)
library(ggplot2)
library(ggsci)
library(stringr)

load('/home/pjy/BRCA/scRNA/Rpro/scRNA/3TF/proj/tnbc.Rdata')
table(Idents(tnbc))

# celltype
pal = pal_d3('category20b')(20)
pal2 = pal_d3('category20c')(20)
DimPlot(tnbc, reduction = "umap", label = T, label.box = T,
        raster = T, pt.size = 0.5,cols = c(pal,pal2)) + NoLegend()
ggsave(filename = '21typecell.pdf', width = 8, height = 7)


# marker
genes_to_check = list(
  `Cancer Cell` = c("WFDC2", "SAA1", "SCGB2A2","C6orf15"),
  `Endothelial Cell` = c("PLVAP","ACKR1","RAMP2","VWF"),
  `Epithelial Cell` = c("CLDN3", "MUCL1", "TFF3","CD24"),
  `Exhausted CD8+ T cell` = c("CD8A" , "CD8B" , "CXCR6"),
  `Fiborblast` = c("DCN","APOD","LUM","SFRP2"),
  `Granulosa Cell` = c("SPP1","RPS20", "FABP4"),
  `Leydig Cell` = c("ASPH","SERPINE1","PFN2"),
  `M1 Macrophage` = c('CCL5','FN1','CXCL9','CCL4','CCL2','CXCL10'),
  `M2 Macrophage` = c("SEPP1","F13A1","FOLR2","SLC40A1","STAB1",'CD163'),
  `Memory T cell` = c("TRAC" , "TRBC2"),
  `Mesenchymal Cell` = c("RGS5","NDUFA4L2"),
  `MKI67+ Progenitor Cell` = c("MKI67","RRM2","UBE2C",'TYMS'),
  `Naive B cell` = c("MS4A1", "CD79B" , "BANK1"),
  `NK cell` = c("GNLY" , "TRDC" ,"NKG7"),
  `NKT Cell` = c("RGCC" , "DNAJB1" ,"DNAJA1"),
  `Neutrophil` = c("S100A8","S100P","PI3","CXCL8"),
  `Plasma cell` = c("IGKV3-15", "IGKC" ,  "IGHG1"),
  `Plasmacytoid dendritic cell` = c('LILRA4',"PLD4",  "NPC2"),
  `Treg cell` = c("BATF",  "TNFRSF18", "TNFRSF1B"),
  `T helper cell` = c("IL7R" ,  "CCR7" ,"CD40LG"),
  `TAM` = c("FCN1","VCAN","EREG","LILRA5","AREG",'CCR2')
)

genes_to_check = lapply(genes_to_check , str_to_upper) # str_to_upper：Converts the case of the string


DotPlot(tnbc,
        features = genes_to_check,
        cols = c("#31A354FF", "#756BB1FF"),
        scale = T,assay='RNA' ) +
  theme(axis.text.x=element_text(angle=45,hjust = 1))
ggsave(filename = '21typecell_marker.pdf', width = 24, height = 7)


# 16gene Vlnplot + FeaturePlot
gene = read.delim(file = '/home/pjy/dock/MINA/mol2/genelist.txt',header = F)
gene = gene$V1

for (i in 1:length(gene)) {
  VlnPlot(tnbc,
          features = gene[i],
          raster = T,
          #pt.size = 0.01,
          cols = c(pal,pal2)) + NoLegend()
  ggsave(filename = paste0('Vlnplot_',gene[i],'.pdf'), width = 7,height = 7)
  
  FeaturePlot(tnbc,
              features = gene[i],
              raster = T,
              #cols = c("lightgrey", "#ff0000", "#00ff00")
  ) 
  ggsave(filename = paste0('FeaturePlot_',gene[i],'.pdf'), width = 7,height = 6)
  
}