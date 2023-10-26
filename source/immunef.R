



suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(reshape2))


Imm1316f <- function(indat,inSet,ingroup,inclinc,insep){
  
  suppressMessages(library(GSVA))
  suppressMessages(library(ComplexHeatmap))
  suppressMessages(library(circlize))
  suppressMessages(library(ggsci))
  suppressMessages(library(GSEABase))
  immdat <- as.matrix(indat)
  geneSets <- inSet
  ingroup <- ingroup
  immclinc <- inclinc
  ssgsea.res <-gsva(expr = immdat,gset.idx.list = geneSets,
                    method = "ssgsea",kcdf = "Gaussian",
                    abs.ranking = T)
  immclinc <- immclinc[order(immclinc[,ingroup]),]
  ssgsea.res <- ssgsea.res[,immclinc$submitter]
  annotation_col = data.frame(Group = immclinc[,ingroup])
  rownames(annotation_col) = immclinc$submitter
  # ann_colors = list(Group = c('low'='#99CC33','high'='#FF9900'))
  Groupscol <- c('#00CCCC','#FF9999')
  names(Groupscol) <- names(table(annotation_col$Group))
  ann_colors = list(Groups = Groupscol)
  mian_col <- colorRampPalette(c("#3300CC", "#3399FF", "white", "#FF3333", "#CC0000"))(50)
  # plotim <-   pheatmap(ssgsea.res, name = "Immune",
  #                      cluster_cols = F,cluster_rows = T,
  #                      annotation_col = annotation_col,
  #                      annotation_colors = ann_colors,
  #                      show_colnames = F,show_rownames = T,
  #                      scale="row",col = mian_col)
  # pdf(file = paste0('fig/step_',insep,'_imm_heatmap.pdf'),width = 8,height = 6)
  # print(plotim)
  # dev.off()
  
  GSVAdat <- as.data.frame(ssgsea.res)
  boxdat16 <- GSVAdat[c(1,4,6,9,11,13,14,16,17,19,22:27),]
  boxdat16 <- rbind(boxdat16,immclinc[,ingroup])
  boxdat16 <- as.data.frame(t(boxdat16))
  names(boxdat16)[17] <- 'RiskGroup'
  boxdat13 <- GSVAdat[-c(1,4,6,9,11,13,14,16,17,19,22:27),]
  boxdat13 <- rbind(boxdat13,immclinc[,ingroup])
  boxdat13 <- as.data.frame(t(boxdat13))
  names(boxdat13)[14] <- 'RiskGroup'
  
  for (immnepy in c(13,16)) {
    work_df <- melt(get(paste0('boxdat',immnepy)), id.vars = "RiskGroup")
    work_df$value <- as.numeric(work_df$value)
    compare_means(value ~ RiskGroup, data = work_df, group.by = "variable")
    names(work_df) <- c('RiskGroup','variable','Score')
    color1 = c('red','green')
    boxpplot <- ggboxplot(work_df, x = "variable", y = "Score",
                          color = "RiskGroup", palette = color1) + 
      stat_compare_means(aes(group = RiskGroup), method = "wilcox.test",label = "p.signif") +
      theme(text = element_text(),
            axis.text.x = element_text(angle = 40,hjust = 1,vjust = 0.98),
            legend.text = element_text())
    ggsave(plot = boxpplot,filename = paste0("fig/step_",insep,"_imm_",immnepy,".pdf"),
           width = 8,height = 5)
  }
  print('successful')
  save(ssgsea.res,file = paste0('rdata/step_',insep,'_ssGSEA.Rdata'))
}



ImmssGSEAf <- function(indat,inSet){
  suppressMessages(library(GSVA))
  suppressMessages(library(GSEABase))
  
  immdat <- indat
  geneSets <- inSet
  # immclinc <- inclinc
  ssgsea.res <-gsva(as.matrix(immdat),geneSets,
                    method = "ssgsea",kcdf = "Gaussian",
                    abs.ranking = T)
  ssgsea.res <- as.data.frame(ssgsea.res)
  
  return(ssgsea.res)
}



Immestif <- function(datpath,inplat,insep){
  suppressMessages(library(estimate))
  
  indar <- datpath
  filterCommonGenes(input.f  = indar,
                    output.f = paste0('tmp/step_',insep,'_estimate.gct'),
                    id       = "GeneSymbol")
  estimateScore(input.ds  = paste0("tmp/step_",insep,"_estimate.gct"), 
                output.ds = paste0("tmp/step_",insep,"_estimate_score.gct"), 
                platform  = inplat)
  scores <- read.table(paste0("tmp/step_",insep,"_estimate_score.gct"),skip = 2,header = T)
  rownames(scores) <- scores[,1]
  scores <- as.data.frame(t(scores))
  scores <- scores[-c(1:2),]
  sampleid <- rownames(scores)
  
  scores <- apply(as.matrix(scores),2,function(x) as.numeric(x))
  scores <- as.data.frame(scores)
  rownames(scores) <- sampleid
  scores$StroGroup <- ifelse(scores$StromalScore > median(scores$StromalScore),'high','low')
  scores$ImmuGroup <- ifelse(scores$ImmuneScore > median(scores$ImmuneScore),'high','low')
  scores$ESTIGroup <- ifelse(scores$ESTIMATEScore > median(scores$ESTIMATEScore),'high','low')
  scores$submitter <- rownames(scores)
  
  return(scores)
}



ImmCibersortf <- function(indat,ingroup,lmdir,insep){
  
  source(paste0(lmdir,'/Cibersort.R'))
  LM22.file <- paste0(lmdir,"/LM22.txt")
  
  # TCGA_exp.file = 'output/step_1_DAT_esti.txt'          
  # load("tmp/step_2_ESTIMATE_SCORE.Rdata")               
  # load("tmp/step_7_TCGA_ROC.Rdata")                     
  # ROC_dat$risk = ifelse(ROC_dat$score > median(ROC_dat$score),'high','low')   
  # 
  # # Gene symbol
  # TCGA_TME.results <- CIBERSORT(LM22.file ,TCGA_exp.file, perm = 50, QN = F)  
  # # perm=1000
  # # 
  # TCGA_TME.results <- TCGA_TME.results[intersect(rownames(TCGA_TME.results),ROC_dat$submitter_id),]
  # # write.csv(TCGA_TME.results, "./output/step_11_CIBERSORT.csv")
  # 
  # ## 2. 
  # # TCGA
  # # group_list <- ifelse(as.numeric(substring(rownames(TCGA_TME.results),14,15)) < 10,
  # #                    "Tumor","Normal") %>% 
  # #  factor(.,levels = c("Normal","Tumor"))
  # 
  # group_list <- ROC_dat$risk %>% factor(.,levels = c("low","high"))
  # table(group_list)
  # 
  # 
  # ## 3. 
  # # 3.1 
  # TME_data <- as.data.frame(TCGA_TME.results[,1:22])
  # TME_data$group <- group_list
  # TME_data$sample <- row.names(TME_data)
  # # 2.2 
  # TME_New = melt(TME_data)
  # ## Using group, sample as id variables
  # colnames(TME_New)=c("Group","Sample","Celltype","Composition")  
  # head(TME_New)
  # # 3.3 
  # plot_order = TME_New[TME_New$Group=="high",] %>% 
  #   group_by(Celltype) %>% 
  #   summarise(m = median(Composition)) %>% 
  #   arrange(desc(m)) %>% 
  #   pull(Celltype)
  # ## `summarise()` ungrouping output (override with `.groups` argument)
  # TME_New$Celltype = factor(TME_New$Celltype,levels = plot_order)
  # box_TME <- ggplot(TME_New, aes(x = Celltype, y = Composition))+ 
  #   labs(y="Cell composition",x= NULL,title = "TME Cell composition")+  
  #   geom_boxplot(aes(fill = Group),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+ 
  #   scale_fill_manual(values = c("#1CB4B8", "#EB7369"))+
  #   theme_classic() + 
  #   stat_compare_means(aes(group =  Group),label = "p.signif",
  #                      method = "wilcox.test",hide.ns = T) + 
  #   theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
  #         axis.title = element_text(size = 12,color ="black"), 
  #         axis.text = element_text(size= 12,color = "black"),
  #         panel.grid.minor.y = element_blank(),
  #         panel.grid.minor.x = element_blank(),
  #         axis.text.x = element_text(angle = 45, hjust = 1 ),
  #         panel.grid=element_blank(),
  #         legend.position = "top",
  #         legend.text = element_text(size= 12),
  #         legend.title= element_text(size= 12))
  # box_TME
  
  ##
  # # 4.1  
  # TCGA_TME_four = as.data.frame(TCGA_TME.results[,1:20])
  # head(TCGA_TME_four,3)
  # # 4.2 
  # # immCell_four_type <- read.table("./database/Cibersort_four_types.txt", header = T, row.names = NULL, sep = "\t")
  # # colnames(TCGA_TME_four) == immCell_four_type$Immune.cells #T
  # # head(immCell_four_type)
  # # 4.3 
  # TCGA_TME_four$group = group_list
  # TCGA_TME_four$sample <- row.names(TCGA_TME_four)
  # TME_four_new = melt(TCGA_TME_four)
  # ## Using group, sample as id variables
  # colnames(TME_four_new) = c("Group","Sample","Immune.cells","Composition")
  # TCGA_TME_four_new2 = left_join(TME_four_new, immCell_four_type, by = "Immune.cells") %>% 
  #   group_by(Sample,Group,Types) %>%
  #   summarize(Sum = sum(Composition))
  # ## `summarise()` regrouping output by 'Sample', 'Group' (override with `.groups` argument)
  # # 
  # box_four_immtypes <- ggplot(TCGA_TME_four_new2, aes(x = Group, y = Sum))+ 
  #   labs(y="Cell composition",x= NULL,title = "TCGA")+  
  #   geom_boxplot(aes(fill = Group),position=position_dodge(0.5),width=0.5,size=0.4,
  #                outlier.alpha = 1, outlier.size = 0.5)+ 
  #   theme_bw() + mytheme + 
  #   scale_fill_manual(values = c("#1CB4B8","#EB7369"))+ 
  #   scale_y_continuous(labels = scales::percent)+
  #   facet_wrap(~ Types,scales = "free",ncol = 4) + 
  #   stat_compare_means(aes(group =  Group),
  #                      label = "p.format",
  #                      method = "wilcox.test",
  #                      size = 3.5,
  #                      hide.ns = T)
  # box_four_immtypes;ggsave("./Output/TCHA_HNSC_Cibersort_four_immune_cell_types.pdf",
  #                          box_four_immtypes ,height= 10,width=25,unit="cm")
}