# DEG from FPKM 

library(tidyverse)
library(limma)
options(stringsAsFactors = F)


# some useful function ----------------------------------------------------
# https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
countToTpm <- function(counts, effLen)
{
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}

countToFpkm <- function(counts, effLen)
{
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}

fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

countToEffCounts <- function(counts, len, effLen)
{
  counts * (len / effLen)
}

tpms <- apply(expMatrix,2,fpkmToTpm)

# lets begin --------------------------------------------------------------
# https://support.bioconductor.org/p/123155/
# https://support.bioconductor.org/p/56275/

gset <- readxl::read_excel('/Users/congliu/OneDrive/kintor/qingiqng/GSE77597_series_matrix.xlsx',
                            skip = 70
                            )

exprSet <- gset %>% 
  filter(gene != '-') %>%
  select(-c(`test_id`,`gene`,locus)) %>% 
  column_to_rownames('gene_id')
Hmisc::describe(exprSet)
exprSet <- exprSet[apply(exprSet,1,function(x){sum(floor(x)==0)<5}),]

# 看看重复的基因
foo <- (exprSet %>% rownames_to_column()) %>% left_join((geset %>% select(gene_id, gene)), 
                                                        by = c('rowname'='gene_id'))
length(sort(table(foo$gene)[table(foo$gene)>1]))

tpms <- apply(exprSet,2,fpkmToTpm)

# limma归一化, 两种归一化的方法需要同时么？
tpms <- log2(tpms + 1)
tpms <- normalizeBetweenArrays(tpms)
boxplot(tpms,las=2)


group_list <- c(rep('Con',2),rep('CA',2),rep('ABA',2), rep('ABTAA',2))

group <- factor(group_list)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
design

#指定那类样本比上那类样本，特别注意有顺序，横杠前的样本比上横杠后的样本
contrast.matrix <- makeContrasts(CA - Con, 
                                 ABA - Con,
                                 ABTAA - Con, 
                                 levels=design)

contrast.matrix

# tpms <- arrayWeights(tpms)
fit <- lmFit(tpms, design)
# fit <- eBayes(fit, trend = TRUE)
# fit <- treat(fit, lfc=log2(1.2), trend=TRUE)
# topTable(fit = fit, coef = 1) # 官网推荐

fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2, trend = TRUE, robust=TRUE)

allDiff1=topTable(fit2,adjust='fdr',coef=1,number=Inf) %>% 
  filter(`adj.P.Val`<0.05)
allDiff2=topTable(fit2,adjust='fdr',coef=2,number=Inf) %>% 
  filter(`adj.P.Val`<0.05)
allDiff3=topTable(fit2,adjust='fdr',coef=3,number=Inf)%>% 
  filter(`adj.P.Val`<0.05)

write_delim(allDiff1,'path/to/save', delim = '\t')





# pca
df=as.data.frame(t(tpms))
df$group=group_list

pp <- autoplot(prcomp( df[,1:(ncol(df)-1)]), 
               data=df,
               colour = 'group')+
  theme_bw()

library("FactoMineR")
library("factoextra") 
dat.pca <- PCA(df[,-ncol(df)], graph = FALSE)

fviz_pca_ind(dat.pca,
             geom.ind = "point",
             col.ind = df$group,
             # palette = c("#00AFBB", "#E7B800"),
             addEllipses = FALSE,
             legend.title = "Groups"
)


# Get the gene lengths and library sizes used to compute the FPKM  --------

# 首先找到参考基因组对应的gene length





