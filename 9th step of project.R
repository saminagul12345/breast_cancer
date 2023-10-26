
source('source/rocf.R')
source('source/BaseTabf.R')
source('source/Staticalf.R')


## 1.load data
load("rdata/step_1_TCGA_EXP.Rdata")
load("rdata/step_5_TCGA_clinc.Rdata")
load("rdata/step_4_DEG_INFO.Rdata")


## 2.Organize data
DEG_INFO <- DEG_INFO[DEG_INFO$State != 'not',]
COX_dat <- TCGA_EXP
COX_clinc <- TCGA_clinc
COX_dat <- COX_dat[rownames(COX_dat) %in% DEG_INFO$SYMBOL,]

COX_tdat <- as.data.frame(t(COX_dat))
COX_tdat$submitter <- rownames(COX_tdat)
TCGA_clinc <- inner_join(COX_clinc,COX_tdat,by = "submitter")
rownames(TCGA_clinc) <- TCGA_clinc$submitter
colnames(TCGA_clinc) <- gsub(colnames(TCGA_clinc), pattern = '[-]', replacement = '_')
rownames(COX_dat) <- gsub(rownames(COX_dat), pattern = '[-]', replacement = '_')


##  Uni variate cox
uncox <- uncoxf(indat = TCGA_clinc,ingene = rownames(COX_dat))
uncox <- uncox[uncox$pvalue < 0.05,]
uncox

## lasso regression
lasso <- lassof(indat = TCGA_clinc,ingene = rownames(uncox),insep = 8)
lasso



# save
write.csv(uncox,file = 'output/step_8_cox_uncox.csv',quote = F)
write.csv(lasso,file = 'output/step_8_cox_lasso.csv',quote = F)
save(lasso,file = 'rdata/step_8_lasso.Rdata')

## train the lasso on merge dataset of metaberic and tcga tnbc

load('rdata/step_1_TCGA_EXP.Rdata')
load('metaberic5_TNBC_exp.Rdata')

tcga<- TCGA_clinc4
metaberic <- metaberic5

merg_mt3 <- merge(tcga,metaberic, all = TRUE)

train <- sample(nrow(merg_mt3),2/3* nrow(merg_mt3))
test = merg_mt3[train]


## Calculate the risk score
names(train_phe)
colnames(train_phe) <- gsub(colnames(train_phe), pattern = '-',replacement = '_')
lasso$gene <- gsub(lasso$gene, pattern = '-', replacement = '_')
train_phe <- theriskf(indat = train_phe,ingene = lasso$gene,incoeff = lasso$coef)
median(train_phe$RiskScore)
table(train_phe$RiskGroup)
save(train_phe,file = 'rdata/step_8_train_phe.Rdata')

## ROC drawing
source('source/rocf.R')
rocf(indat = train_phe)

##  High and low risk groups KM
sfit <- survfit(Surv(OS.time,OS) ~ RiskGroup, data = train_phe)
plot4 <- ggsurvplot(sfit,palette = c("#ea0102","#63b931"),
                    risk.table = F,
                    pval =TRUE,pval.size= 4.3,pval.coord=c(1,0.03),
                    pval.method = F,
                    xlim = c(0,150),
                    conf.int = T,xlab ="Time (months)",ylab = 'Survival probability',
                    ggtheme = theme_test(),
                    legend = c(0.82,0.83),
                    legend.title = "Risk Group",legend.labs = c("High risk", "Low risk"))
cox<-coxph(Surv(OS.time,OS) ~ RiskScore,data = train_phe)
res<-summary(cox)
plot2 <-   plot4$plot +
  annotate("text",x=11,y=0.24,label=paste("HR = ",signif(res$conf.int[1],2))) +
  annotate("text",x=20,y=0.17,label=paste("95%CI = ",round(res$conf.int[3],2),"-",
                                          signif(res$conf.int[4],2)))+
  annotate("text",x=20,y=0.1,label=paste("Cox P = ",signif(res$coefficients[5],3)))
print(plot2)
pdf(file = 'fig/step_90_Mdule_KMinterST.pdf',width = 4,height = 4)
print(plot2)
dev.off()


##  Risk factor correlation chart
source('source/factorf.R')
factorf(indat = train_phe,inGene = lasso$gene,insep = '8_TCGA10')