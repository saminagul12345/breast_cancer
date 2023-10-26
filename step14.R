rm(list = ls())
options(stringsAsFactors = F)
gc()

library(Seurat)
library(ggplot2)
load('/home/pjy/BRCA/scRNA/Rpro/scRNA/3TF/proj/tnbc.Rdata')
table(Idents(tnbc))

Epi = tnbc[,Idents(tnbc) %in% 'Epithelial Cell']
Epi

Epi <- NormalizeData(Epi)
Epi <- FindVariableFeatures(Epi, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Epi)
Epi <- ScaleData(Epi, features = all.genes)
Epi <- RunPCA(Epi, features = VariableFeatures(object = Epi))
Epi <- JackStraw(Epi, num.replicate = 100) 
Epi <- ScoreJackStraw(Epi, dims = 1:20)
JackStrawPlot(Epi, dims = 1:20)
ggsave(filename = "Epi_JackPlot.png",width = 16,height = 10,path = "/home/datahup/mina/scRNA/R_data/my_Rscript/geo/fig4/Epi")
ElbowPlot(Epi,ndims = 50) 
ggsave(filename = "Epi_ElbowPlot.png",width = 12,height = 10,path = "/home/datahup/mina/scRNA/R_data/my_Rscript/geo/fig4/Epi")

library(clustree)
Epi <- FindNeighbors(Epi, dims = 1:12)
Epi <- FindClusters(object = Epi,
                        resolution = c(seq(0,1,by = 0.1)))
clustree(Epi@meta.data, prefix = "RNA_snn_res.") 
ggsave(filename = "Epi_resolution(0-1).png",width = 20,height = 14,path = "/home/datahup/mina/scRNA/R_data/my_Rscript/geo/fig4/Epi")
Idents(object = Epi) <- "RNA_snn_res.0.3"
Epi@meta.data$seurat_clusters = Epi@meta.data$RNA_snn_res.0.3
head(Idents(Epi), 5)#查看前5个细胞的分类ID
# UMAP
Epi <- RunUMAP(Epi, dims = 1:12)
DimPlot(Epi, reduction = "umap", label = TRUE,raster=FALSE)
ggsave(filename = "Epi_UMAP-label.pdf",width = 8,height = 6,path = "/home/datahup/mina/scRNA/R_data/my_Rscript/geo/fig4/Epi")
save(Epi,file = 'Epi.Rdata')
load(file = 'Epi.Rdata')
# Step1 Seurat -> monocle
data <- as(as.matrix(Epi@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = Epi@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)

monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
print(head(fData(monocle_cds)))


HSMM=monocle_cds
disp_table <- dispersionTable(HSMM)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
plot_ordering_genes(HSMM)
plot_pc_variance_explained(HSMM, return_all = F)

HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 6,
                        reduction_method = 'tSNE', verbose = T)
HSMM <- clusterCells(HSMM, num_clusters = 10)
plot_cell_clusters(HSMM)

disp_table <- dispersionTable(HSMM)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
HSMM <- setOrderingFilter(HSMM, disp.genes)
plot_ordering_genes(HSMM)

HSMM <- reduceDimension(HSMM, max_components = 3,
                        method = 'DDRTree')
HSMM <- orderCells(HSMM)
plot_cell_trajectory(HSMM, color_by = "seurat_clusters")+
  facet_wrap(~seurat_clusters, nrow = 1)
ggsave(filename = 'seurat_clusters.png',width = 10,height = 5,path = '/home/datahup/mina/scRNA/R_data/my_Rscript/geo/fig4/Epi')


plot_cell_trajectory(HSMM, color_by = "State")

plot_cell_trajectory(HSMM, color_by = "Pseudotime")

plot_cell_trajectory(HSMM, color_by = "State") +
  facet_wrap(~State, nrow = 1)
ggsave(filename = 'State.png',width = 10,height = 5,path = '/home/datahup/mina/scRNA/R_data/my_Rscript/geo/fig4/Epi')

library(dplyr)
Time_diffE <- differentialGeneTest(HSMM[disp.genes,], cores = 4, 
                                   fullModelFormulaStr = "~sm.ns(Pseudotime)")
# Time_diff <- Time_diff[,c(5,2,3,4,1,6,7)] #ægeneæ¾åé¢ï¼ä¹å¯ä»¥ä¸æ¹
write.csv(Time_diffE, "Time_diff_E.csv", row.names = F)
Time_genes <- Time_diffE %>% pull(gene_short_name) %>% as.character()
# num_clustersä¸ºäººä¸ºè®¾ç½®çèç±»
p=plot_pseudotime_heatmap(HSMM[Time_genes,], num_clusters=3, show_rownames=T, return_heatmap=T)

ggsave("Time_heatmapAll.png", p, width = 10, height = 12,path = "/home/datahup/mina/scRNA/R_data/my_Rscript/geo/fig4/Epi")

p$tree_row
#Call:
 # hclust(d = d, method = method)

#Cluster method   : ward.D2 
# Number of objects: 1075
clusters <- cutree(p$tree_row, k = 3)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
table(clustering)
#Gene_Clusters
#1   2   3 
#612 222 241
write.csv(clustering, "Time_clustering_E.csv", row.names = F)
library(ggpubr)
library(Seurat)
df <- pData(HSMM) 
table(df$State)
df$NewState <- ifelse(df$State == '1', 'State1',
                      ifelse(df$State == '2','State2',
                             ifelse(df$State == '3', 'State3',    
                                    ifelse(df$State == '4', 'State3',
                                           ifelse(df$State == '5', 'State4',
                                                  ifelse(df$State == '6', 'State5',
                                                         ifelse(df$State == '7', 'State6',
                                                              
                                                                              ifelse(df$State =='4','state3','State5'))))))))

table(df$NewState)
State <- data.frame(State = df$NewState, row.names = rownames(df), Cell = rownames(df))
State <- State[order(State$State),]

table(State$State)
Epi <- AddMetaData(Epi, metadata = State)
Idents(object = Epi) <- "State"
head(Idents(Epi), 5)#æ¥çå5ä¸ªç»èçåç±»ID
table(Idents(Epi))

Statemarker <- FindAllMarkers(Epi)
table(Statemarker$cluster)
top50 <-Statemarker %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
DoHeatmap(Epi, features = top50$gene) + NoLegend()
ggsave(filename = 'State_Marker50.png',width = 8,height = 12,path = '/home/datahup/mina/scRNA/R_data/my_Rscript/geo/fig4/Epi')
save(Statemarker,file = 'Epi_Statemarker.Rdata')

library(dplyr)
gene0 <- unique(Statemarker$gene)
Time_diff <- differentialGeneTest(HSMM[gene0,], cores = 4, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
# num_clustersä¸ºäººä¸ºè®¾ç½®çèç±»
p <- plot_pseudotime_heatmap(HSMM[Time_genes,], num_clusters = 4,show_rownames=T, return_heatmap=T)
ggsave(filename = 'State_Marker_pseudotime.png',plot = p,width = 8,height = 12,path = '/home/datahup/mina/scRNA/R_data/my_Rscript/geo/fig4/Epi')

HSMM@phenoData@data$State <- df$NewState
p1 <- plot_cell_trajectory(HSMM, color_by = "seurat_clusters")
p2 <- plot_cell_trajectory(HSMM, color_by = "State", show_branch_points = FALSE)
p3 <- plot_cell_trajectory(HSMM, color_by = "Pseudotime",show_branch_points = FALSE)
plotc <- p1|p2|p3
plotc
ggsave(filename = 'Clu+State+Tra.png',width = 12,height = 5,plot = plotc,path = '/home/datahup/mina/scRNA/R_data/my_Rscript/geo/fig4/Epi')


plot_cell_trajectory(HSMM,color_by = "State")+facet_wrap(~State,nrow=1)
ggsave(filename = 'State.png',width = 10,height = 5,path = '/home/datahup/mina/scRNA/R_data/my_Rscript/geo/fig4/Epi')
my_genes <- row.names(subset(fData(HSMM),
                             gene_short_name %in% c('BMP4', 'CCBE1', 'CELSR3', 'CT83', 'CXCL11', 'EGR2', 'GLDC', 'GPRC5C', 
                                                    'TRO', 'STMN2', 'SCGB2A2', 'RUNDC3B', 'PROS1', 'PCDHGA3', 'IL1RL1', 'UGT2B1')))
cds_subset <- HSMM[my_genes,]


plot_genes_in_pseudotime(cds_subset, color_by =  "State")
plot_genes_jitter(cds_subset, grouping = "State", color_by = "State")
plot_genes_violin(cds_subset, grouping = "State", color_by = "State")

save(Epi,HSMM,file = "mono_Epi.Rdata")


library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
GO_GSEA_list <- list()
KEGG_GSEA_list <- list()

for (i in 1:4) {
  
  name <- paste0('State',i)
  dat <- subset(Statemarker,Statemarker$cluster == name)
  dat$SYMBOL <- dat$gene
  
  ENTREZID1 <- bitr(unique(dat$SYMBOL), fromType = "SYMBOL",
                    toType = c("ENTREZID"),
                    OrgDb = org.Hs.eg.db)
  
  dat <- subset(dat,dat$SYMBOL %in% ENTREZID1$SYMBOL)
  dat <- merge(dat,ENTREZID1)
  geneList=dat$avg_log2FC
  names(geneList)=dat$ENTREZID
  geneList=sort(geneList,decreasing = T)
  
  #GO GSEA
  GO <- gseGO(
    geneList, #gene_fc
    ont = "ALL",# "BP"ã"MF"å"CC"æ"ALL"
    OrgDb = org.Hs.eg.db,#äººç±»æ³¨éåºå 
    keyType = "ENTREZID",
    pvalueCutoff = 0.9,
    pAdjustMethod = "BH")#på¼æ ¡æ­£æ¹æ³
  
  sortGO<-GO[order(GO$enrichmentScore, decreasing = T),]#æç§enrichment scoreä»é«å°ä½æåº
  head(sortGO)
  dim(sortGO)
  # write.table(sortGO,paste0(paste0('State',i),'_','gsea_sortGO.txt')) #ä¿å­ç»æ
  GO_GSEA_list[[i]] <- sortGO
  GO_GSEA_list[[i+4]] <- GO
  
  #KEGG GSEA
  KEGG <- gseKEGG(
    geneList,
    organism = "hsa",
    keyType = "kegg",
    pvalueCutoff = 0.9,
    pAdjustMethod = "BH")
  
  
  sortKEGG<-KEGG[order(KEGG$enrichmentScore, decreasing = T),]#æç§enrichment scoreä»é«å°ä½æåº
  head(sortKEGG)
  dim(sortKEGG)
  # write.table(sortKEGG,paste0(paste0('State',i),'_','gsea_sortKEGG.txt')) #ä¿å­ç»æ
  KEGG_GSEA_list[[i]] <- sortKEGG
  KEGG_GSEA_list[[i+4]] <- KEGG
  
}

save(KEGG_GSEA_list,GO_GSEA_list,file = 'State_GSEA_E.Rdata')
load(file = 'State_GSEA_E.Rdata')

library(enrichplot)
library(dplyr)
library(ggplot2)
library(ggsci)
pal = pal_ucscgb()(20)
State1_KEGG <- as.data.frame(KEGG_GSEA_list[1]) %>% filter(pvalue < 0.05)
write.csv(State1_KEGG,file = 'State1E_KEGG.csv')
paths <- State1_KEGG$ID#éåä½ éè¦å±ç¤ºçéè·¯ID
gseaplot2(KEGG_GSEA_list[[5]],paths, pvalue_table = TRUE,title = 'State1 KEGG Enrichment',
          rel_heights = c(2, 0.3,0.3),subplots = c(1:3), color = pal[1:9])
ggsave(filename = 'State1 KEGG Enrichment.png',width = 16,height = 20,path = '/home/datahup/mina/scRNA/R_data/my_Rscript/geo/fig4/Epi')


State2_KEGG <- as.data.frame(KEGG_GSEA_list[2]) %>% filter(pvalue < 0.05)
write.csv(State2_KEGG,file = 'State2E_KEGG.csv')
paths <- State2_KEGG$ID#éåä½ éè¦å±ç¤ºçéè·¯ID
gseaplot2(KEGG_GSEA_list[[6]],paths, pvalue_table = TRUE,title = 'State2 KEGG Enrichment',
          rel_heights = c(2, 0.3,0.3),subplots = c(1:3), color = pal[1:2])
# paths <- c("hsa05171", "hsa03010")#éåä½ éè¦å±ç¤ºçéè·¯ID
# gseaplot2(KEGG_GSEA_list[[6]],paths, pvalue_table = TRUE,title = 'State2 KEGG Enrichment')
ggsave(filename = 'State2 KEGG Enrichment.png',width = 16,height = 20,path = '/home/datahup/mina/scRNA/R_data/my_Rscript/geo/fig4/Epi')


State3_KEGG <- as.data.frame(KEGG_GSEA_list[3]) %>% filter(pvalue < 0.05)
write.csv(State3_KEGG,file = 'State3E_KEGG.csv')
paths <- State3_KEGG$ID#éåä½ éè¦å±ç¤ºçéè·¯ID
gseaplot2(KEGG_GSEA_list[[7]],paths, pvalue_table = TRUE,title = 'State3 KEGG Enrichment',
          rel_heights = c(2, 0.3,0.3),subplots = c(1:3), color = pal[1:8])
# paths <- c("hsa05140", "hsa04613","hsa03010","hsa05171")#éåä½ éè¦å±ç¤ºçéè·¯ID
# gseaplot2(KEGG_GSEA_list[[7]],paths, pvalue_table = TRUE,title = 'State3 KEGG Enrichment')
ggsave(filename = 'State3 KEGG Enrichment.png',width = 18,height = 20,path = '/home/datahup/mina/scRNA/R_data/my_Rscript/geo/fig4/Epi')


State4_KEGG <- as.data.frame(KEGG_GSEA_list[4]) %>% filter(pvalue < 0.05)
write.csv(State4_KEGG,file = 'State4E_KEGG.csv')
paths <- State4_KEGG$ID#éåä½ éè¦å±ç¤ºçéè·¯ID
gseaplot2(KEGG_GSEA_list[[8]],paths, pvalue_table = TRUE,title = 'State4 KEGG Enrichment',
          rel_heights = c(2, 0.3,0.3),subplots = c(1:3), color = pal[1:6])
# paths <- c("hsa05323", "hsa05171","hsa04144","hsa03010")#éåä½ éè¦å±ç¤ºçéè·¯ID
# gseaplot2(KEGG_GSEA_list[[8]],paths, pvalue_table = TRUE,title = 'State4 KEGG Enrichment')
ggsave(filename = 'State4 KEGG Enrichment.png',width = 16,height = 20,path = '/home/datahup/mina/scRNA/R_data/my_Rscript/geo/fig4/Epi')

State5_KEGG <- as.data.frame(KEGG_GSEA_list[4]) %>% filter(pvalue < 0.05)
write.csv(State5_KEGG,file = 'State5E_KEGG.csv')
paths <- State5_KEGG$ID#éåä½ éè¦å±ç¤ºçéè·¯ID
gseaplot2(KEGG_GSEA_list[[8]],paths, pvalue_table = TRUE,title = 'State5 KEGG Enrichment',
          rel_heights = c(2, 0.3,0.3),subplots = c(1:3), color = pal[1:6])
# paths <- c("hsa05323", "hsa05171","hsa04144","hsa03010")#éåä½ éè¦å±ç¤ºçéè·¯ID
# gseaplot2(KEGG_GSEA_list[[8]],paths, pvalue_table = TRUE,title = 'State4 KEGG Enrichment')
ggsave(filename = 'State5 KEGG Enrichment.png',width = 16,height = 20,path = '/home/datahup/mina/scRNA/R_data/my_Rscript/geo/fig4/Epi')

State6_KEGG <- as.data.frame(KEGG_GSEA_list[4]) %>% filter(pvalue < 0.05)
write.csv(State6_KEGG,file = 'State6E_KEGG.csv')
paths <- State6_KEGG$ID#éåä½ éè¦å±ç¤ºçéè·¯ID
gseaplot2(KEGG_GSEA_list[[8]],paths, pvalue_table = TRUE,title = 'State6 KEGG Enrichment',
          rel_heights = c(2, 0.3,0.3),subplots = c(1:3), color = pal[1:6])
 paths <- c("hsa01200", "hsa00270","hsa00630","hsa00620")#éåä½ éè¦å±ç¤ºçéè·¯ID
gseaplot2(KEGG_GSEA_list[[8]],paths, pvalue_table = TRUE,title = 'State6 KEGG Enrichment')
ggsave(filename = 'State6 KEGG Enrichment.png',width = 16,height = 20,path = '/home/datahup/mina/scRNA/R_data/my_Rscript/geo/fig4/Epi')




load(file = 'mono_Epi.Rdata')

celltype = as.data.frame(table(Idents(Epi)))
celltype1 = as.character(celltype$Var1)
all_marker = FindAllMarkers(Epi, only.pos = TRUE, min.pct = 0.35, logfc.threshold = 0.5)
all_marker.top10 = all_marker %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
table(duplicated(all_marker.top10$gene))
cell_dat = list()
for (i in celltype1) {
  cell = Epi[,Idents(Epi) %in% i]
  dat = as.matrix(cell@assays[["RNA"]]@data) # 提取基因表达数据
  dat = dat[all_marker.top10$gene,]
  dat1 = data.frame(gene = rownames(dat), 
                    exp = rowSums(dat)) # 计算该簇中所有基因的总和
  names(dat1)[2] = i
  cell_dat[[i]] = dat1
}
cell_dat1 = do.call(cbind, cell_dat)
cell_dat2 = cell_dat1[,-c(3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41)]
names(cell_dat2) = c('Gene', 'State3', 'State5', 'State4', 'State2',
                     'State6', 'State1')
table(duplicated(cell_dat2$Gene))
cell_dat2 = subset(cell_dat2, !duplicated(cell_dat2$Gene))
write.table(cell_dat2,file="6state_top10.txt",sep = "\t")# EXCEL打开，删除第一列，第一行右移一格
 X <- read.table(file = '6state.txt',header=T,sep="\t",row.names=1,check.names=F)  # 检查数据

table(phe2$SubType)
phe3 = subset(phe2, phe2$SubType == 'Basal' )
phe4 = subset(phe2, phe2$SubType == 'Normal' )
phe5 = rbind(phe3,phe4)
Basal_dat = dat[,rownames(phe5)]
Basal_dat[1:6,1:6]
table(duplicated(rownames(Basal_dat)))
Basal_dat=as.matrix(Basal_dat)
dimnames=list(rownames(Basal_dat),colnames(Basal_dat))
data=matrix(as.numeric(as.matrix(Basal_dat)),nrow=nrow(Basal_dat),dimnames=dimnames)
library(limma)
data=avereps(data)
data=data[rowMeans(data)>0,]
data[data<0] = NA
data = na.omit(data)
is.recursive(data)
is.atomic(data)
dim(data)
write.csv(data,file = "CIBERSORT_Basal_guolv.csv")
a = read.delim(file = "CIBERSORT_Basal_guolv.txt")
a[1:6,1:6]
table(duplicated(a[,1]))
a = subset(a, !duplicated(a[,1]))
write.table(a,file = "CIBERSORT_Basal_guolv.txt",sep = '\t')
#X <- read.table(file = 'CIBERSORT_Basal_guolv.txt',header=T,sep="\t",row.names=1,check.names=F)  # 检查数据
# X[1:4,1:4]
#                 TCGA.A1.A0SK.01A TCGA.A1.A0SO.01A TCGA.A1.A0SP.01A TCGA.A2.A04P.01A
# DDX11L1              0.000000         1.000000         1.000000         1.584963
# WASH7P               5.169925         4.247928         5.392317         6.554589
# MIR6859-3            0.000000         0.000000         0.000000         1.000000
# RP11-34P13.3         0.000000         0.000000         0.000000         0.000000

source("Cibersort.R")
library(preprocessCore)
library(e1071)
library(parallel)
result1 <- CIBERSORT("6state.txt", "CIBERSORT_Basal_guolv.txt", perm = 100, QN = T)
save(result1,file = 'CIBERSORT_top10_result.Rdata')

load(file = 'CIBERSORT_top10_result.Rdata')
re <- result1[,-(7:9)] # tnbc+normal
re1 = subset(re, rownames(re) %in% rownames(phe3))

library(RColorBrewer)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
dat <- re %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Cell_type,value = Proportion,-Sample)
dat1 <- re1 %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Cell_type,value = Proportion,-Sample)
library(ggsci)
mypalette = pal_jco()(10)
mypalette1 = pal_ucscgb()(10)
mypalette2 = pal_lancet()(8)
ggplot(dat1,aes(Sample,Proportion,fill = Cell_type)) + 
  geom_bar(stat = "identity") +
  labs(fill = "Cell Type",x = "",y = "Cibersort Proportion") + 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = c(mypalette1,mypalette,mypalette2))
ggsave(filename = "Cibersort_bar.pdf",width = 12,height = 8, path = '/home/datahup/mina/scRNA/R_data/my_Rscript/geo/fig4/Epi')
a = dat1 %>% 
  group_by(Cell_type) %>% 
  summarise(m = median(Proportion)) %>% 
  arrange(desc(m)) %>% 
  pull(Cell_type)
dat1$Cell_type = factor(dat1$Cell_type,levels = a)
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
ggplot(dat1,aes(Cell_type,Proportion,fill = Cell_type)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Cibersort Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_fill_manual(values = mypalette(28))
ggsave(filename = "Cibersort_box.pdf",width = 11,height = 6, path = '/home/datahup/mina/scRNA/R_data/my_Rscript/geo/fig4/Epi')


top10 <-Statemarker %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Epi, features = top10$gene) + NoLegend()
ggsave(filename = 'State_Marker10.png',width = 12,height = 15,path = '/home/datahup/mina/scRNA/R_data/my_Rscript/geo/fig4/Epi')

