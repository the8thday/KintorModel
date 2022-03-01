# cell miner

# analysis for cell miner website

library(readxl)
library(tidyverse)
library(magrittr)
library(impute)
library(correlation)
options(stringsAsFactors = F)


rt <- read_excel('/Users/congliu/OneDrive/kintor/nci60/DTP_NCI60_ZSCORE.xlsx',
                 skip = 8
                 )
i_cells <- Hmisc::Cs(
  LE:CCRF_CEM,
  LE:HL_60,
  LE:K_562,
  LE:MOLT_4,
  LE:RPMI_8226,
  LE:SR,
  LC:A549,
  LC:EKVX,
  LC:HOP_62,
  LC:HOP_92,
  LC:NCI_H226,
  LC:NCI_H23,
  LC:NCI_H322M,
  LC:NCI_H460,
  LC:NCI_H522
)
i_cells <- gsub('_', '-', i_cells)

rt %<>% filter(`FDA status` %in% c('Clinical trial', 'FDA approved')) %>%
  dplyr::select(-`Total experiments`, -`Total after quality control`)
drug <- rt[,-c(2:6)] %>%
  dplyr::select(`NSC # b`, starts_with(c('LC:','LE:'))) %>%
  # select(any_of(i_cells))
  column_to_rownames('NSC # b')

rna  <-  read_excel('/Users/congliu/OneDrive/kintor/nci60/RNA__RNA_seq_composite_expression.xls',
                    skip = 10
                    )

rna <- rna[,-c(2:6)] %>%
  dplyr::select(`Gene name d`, starts_with(c('LC:','LE:'))) %>%
  column_to_rownames('Gene name d')

sum(is.na(drug))
dimnames <- list(rownames(drug),colnames(drug))
data <- matrix(as.numeric(as.matrix(drug)),nrow=nrow(drug),dimnames=dimnames)
mat <- impute.knn(data = data, colmax = 1)
drug <- limma::avereps(mat$data)


genelist <- c('TP53','MYC','BCL2','BCL6','KRAS','JAK1','JAK2','FGFR1','FGFR2',
              'FGFR3','RUNX1','CREBBP','EP300','MALT1','BCR')
genelist <- intersect(genelist,row.names(rna))
rna <- rna[genelist,]

# 分别计算每个基因表达与不同药物之间的Pearson相关系数
# 据此可以看到和我们相似的药物在这些细胞系中和哪些基因表达相关
outTab <- data.frame()

for(Gene in row.names(rna)){
  x <- as.numeric(rna[Gene,])
  for(Drug in row.names(drug)){
    y <- as.numeric(drug[Drug,])
    corT <- cor.test(x,y,method="pearson")
    cor <- corT$estimate
    pvalue <- corT$p.value
    if(pvalue < 0.01){
      outVector <- cbind(Gene,Drug,cor,pvalue)
      outTab <- rbind(outTab,outVector)
    }
  }
}

outTab <- outTab[order(as.numeric(as.vector(outTab$pvalue))),]
write_excel_csv(outTab, file = '/Users/congliu/OneDrive/kintor/nci60/cor_res.csv')

myc_drug <- outTab %>% filter(Gene == 'MYC') %>% pull(Drug)
# read_excel('/Users/congliu/OneDrive/kintor/nci60/DTP_NCI60_ZSCORE.xlsx',
#            skip = 8
# ) %>% janitor::clean_names() %>% filter(nsc_number_b %in% myc_drug) %>%
#   write_excel_csv(file = '/Users/congliu/OneDrive/kintor/nci60/MYC_drug.csv')


# plot，展示基因表达和药物IC50之间的关系

library(ggpubr)

plotList_1 <- list()
corPlotNum <- 16
if(nrow(outTab)<corPlotNum){
  corPlotNum=nrow(outTab)
}

for(i in 1:corPlotNum){
  Gene <- outTab[i,1]
  Drug <- outTab[i,2]
  x <- as.numeric(rna[Gene,])
  y <- as.numeric(drug[Drug,])
  cor <- sprintf("%.03f",as.numeric(outTab[i,3]))
  pvalue=0
  if(as.numeric(outTab[i,4])<0.001){
    pvalue="p<0.001"
  }else{
    pvalue=paste0("p=",sprintf("%.03f",as.numeric(outTab[i,4])))
  }
  df1 <- as.data.frame(cbind(x,y))
  p1=ggplot(data = df1, aes(x = x, y = y))+
    geom_point(size=1)+
    stat_smooth(method="lm",se=FALSE, formula=y~x)+
    labs(x="Expression",y="IC50",title = paste0(Gene,", ",Drug),subtitle = paste0("Cor=",cor,", ",pvalue))+
    theme(axis.ticks = element_blank(), axis.text.y = element_blank(),axis.text.x = element_blank())+
    theme_bw()
  plotList_1[[i]]=p1
}

plotList_2 <- list()
corPlotNum <- 16
if(nrow(outTab)<corPlotNum){
  corPlotNum=nrow(outTab)
}


for(i in 1:corPlotNum){
  Gene <- outTab[i,1]
  Drug <- outTab[i,2]
  x <- as.numeric(rna[Gene,])
  y <- as.numeric(drug[Drug,])
  df1 <- as.data.frame(cbind(x,y))
  colnames(df1)[2] <- "IC50"
  df1$group <- ifelse(df1$x > median(df1$x), "high", "low")
  compaired <- list(c("low", "high"))
  p1 <- ggboxplot(df1,
                  x = "group", y = "IC50",
                  fill = "group", palette = c("#00AFBB", "#E7B800"),
                  add = "jitter", size = 0.5,
                  xlab = paste0("The_expression_of_", Gene),
                  ylab = paste0("IC50_", Drug)) +
    stat_compare_means(comparisons = compaired,
                       method = "wilcox.test",   #设置统计方法
                       symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "**", "*", "ns")))
  plotList_2[[i]]=p1
}



nrow <- ceiling(sqrt(corPlotNum))
ncol <- ceiling(corPlotNum/nrow)
bar <- ggarrange(plotlist=plotList_1,nrow=nrow,ncol=ncol)
foo <- ggarrange(plotlist=plotList_2,nrow=nrow,ncol=ncol)

ggsave(
  plot = foo,
  filename = '/Users/congliu/OneDrive/kintor/nci60/box.pdf',
  height = 15,
  width = 15
)

ggsave(
  plot = bar,
  filename = '/Users/congliu/OneDrive/kintor/nci60/line.pdf',
  height = 15,
  width = 15
)
