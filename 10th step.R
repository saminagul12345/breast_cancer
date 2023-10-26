source('source/rocf.R')
source('source/oclr.R')
source('source/BaseTabf.R')
source('source/Staticalf.R')
source('source/Id2Symbolf.R')
source('source/factorf.R')


# load data
lasso <- read.csv('output/step_8_cox_lasso.csv',row.names = 1)



# ready COX data
COX_dat <- test_exp
COX_dat <- COX_dat[rownames(COX_dat) %in% lasso$gene,]
COX_clinc <- test_phe

COX_tdat <- as.data.frame(t(COX_dat))
COX_tdat$id <- rownames(COX_tdat)
COX_mdat <- inner_join(COX_clinc,COX_tdat,by = "id")
rownames(COX_mdat) <- COX_mdat$id
COX_mdat$OS.time <- as.numeric(COX_mdat$OS.time)
COX_mdat$OS <- as.numeric(COX_mdat$OS)


## Calculate the score 
Risk_mdat <- theriskf(indat = COX_mdat,ingene = lasso$gene,incoeff = lasso$coef)
Risk_mdat$RiskScore
median(Risk_mdat$RiskScore)
table(Risk_mdat$RiskGroup)

## ROC curve
rocf(indat = Risk_mdat)

## High and low risk groups KM
sfit <- survfit(Surv(OS.time,OS) ~ RiskGroup, data = Risk_mdat)
plot4 <- ggsurvplot(sfit,palette = c("#ea0102","#63b931"),
                    risk.table = F,
                    pval =TRUE,pval.size= 4.3,pval.coord=c(1,0.03),
                    pval.method = F,
                    xlim = c(0,150),
                    conf.int = T,xlab ="Time (months)",ylab = 'Survival probability',
                    ggtheme = theme_test(),
                    legend = c(0.82,0.83),
                    legend.title = "Risk Group",legend.labs = c("High risk", "Low risk"))
cox<-coxph(Surv(OS.time,OS) ~ RiskScore,data = Risk_mdat)
res<-summary(cox)
plot2 <-   plot4$plot +
  annotate("text",x=11,y=0.24,label=paste("HR = ",signif(res$conf.int[1],2))) +
  annotate("text",x=20,y=0.17,label=paste("95%CI = ",round(res$conf.int[3],2),"-",
                                          signif(res$conf.int[4],2)))+
  annotate("text",x=20,y=0.1,label=paste("Cox P = ",signif(res$coefficients[5],3)))
print(plot2)
pdf(file = 'fig/step_90_Mdule_KMinterSTest.pdf',width = 4,height = 4)
print(plot2)
dev.off()