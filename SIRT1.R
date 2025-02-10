setwd("E:\\于恒响2023\\USP7\\生信\\SLE生信\\GSE65391 SLE对")
load(".Rdata")
# 导入所需的R包
library(data.table)
library(ggstatsplot)
library(ggplot2)
library(ggcorrplot)
# 读取表达矩阵文件
exprSet <- read.csv("GSE46239_data_matrix.csv", row.names = 1)
#筛选#
#exprSet <- exprSet[ ,colnames(exprSet)[grep("SC2", colnames(exprSet))]]根据列名提取#
#write.table(exprSet, file = "exprSet.txt", sep = "\t", row.names = T, col.names = F)#
#exprSet <- exprSet[,seq(1,ncol(exprSet),2)]筛选奇数列#
#exprSet <- exprSet[, -grep("control", colnames(exprSet))]筛选包含？列删除#



###两基因相关性##########################
#转置
dat<-as.data.frame(t(exprSet))
# 保留?行，病人
datSLE <- dat[15:69, ]
#后30个样本，正常
datHD <- dat[(nrow(dat) - 29):nrow(dat), ]
#两基因相关性
ggscatterstats(data = dat, 
               y = SIRT1, 
               x = IFNB1,
               centrality.para = "mean",                              
               margins = "both",                                         
               xfill = "#CC79A7", 
               yfill = "#009E73", 
               marginal.type = "histogram",
               title = "GSE46239")
ggsave("SIRT1 IFNB1.tiff", device = tiff, width = 8, height = 6)


###多基因相关性#####  
# 向量sirtGenes包含SIRT基因的名称  
sirtGenes <- c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7")  
# 向量ifnGenes包含IFN基因的名称  
ifnGenes <- c("IFNB1","IFNA5","IFNA4","IFNA2","IFNA1")  
#"IFNA6", "IFNA7", "IFNA8", "IFNA10", "IFNA13", "IFNA14", "IFNA16", "IFNA17", "IFNA21",#

# 创建一个空的相关性矩阵  
# 矩阵的行数和列数分别等于sirtGenes和ifnGenes的长度  
corMatrix <- matrix(NA, nrow = length(sirtGenes), ncol = length(ifnGenes))  

# 计算每个"SIRT"基因与每个"IFN"基因之间的相关性  
# 使用两个嵌套的for循环遍历sirtGenes和ifnGenes中的所有基因  
for (i in 1:length(sirtGenes)) {  
  for (j in 1:length(ifnGenes)) {  
    # 获取当前遍历的SIRT基因名称  
    sirtGene <- sirtGenes[i]  
    # 获取当前遍历的IFN基因名称  
    ifnGene <- ifnGenes[j]  
    # 从exprSet中获取当前SIRT基因的数据  
    sirtData <- dat[,sirtGene]  
    # 从exprSet中获取当前IFN基因的数据  
    ifnData <- dat[,ifnGene]  
    # 计算当前SIRT基因和当前IFN基因之间的相关性，并将结果存入corMatrix矩阵的对应位置  
    corMatrix[i, j] <- cor(sirtData, ifnData)  
  }  
}  

# 绘制相关性矩阵热图  
# 设置corMatrix矩阵的行名和列名为对应的基因名称  
rownames(corMatrix) <- sirtGenes  
colnames(corMatrix) <- ifnGenes  
# 使用corrplot函数绘制相关性矩阵的热图，使用颜色表示相关性的强弱  
ggcorrplot(corMatrix,tl.cex = 12,lab_size = 8)
ggsave("SIRT1-7 IFN.tiff", device = tiff, width = 8, height = 6)

### 保留前?列，病人####################
exprSet <- exprSet[, 1:25]
#批量计算相关性
batch_cor <- function(gene){
  y <- as.numeric(exprSet[gene,])
  rownames <- rownames(exprSet)
  do.call(rbind,future_lapply(rownames, function(x){
    dd  <- cor.test(as.numeric(exprSet[x,]),y,type="spearman")
    data.frame(gene=gene,mRNAs=x,cor=dd$estimate,p.value=dd$p.value )
  }))
}

library(future.apply)
plan(multicore)
system.time(dd <- batch_cor("SIRT1"))

write.table(dd, file = "USP7相关性.txt", sep = "\t", row.names = FALSE)

#GSEA分析#
gene <- dd$mRNAs
## 转换
library(clusterProfiler)
gene = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
## 去重
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)

gene_df <- data.frame(logFC=dd$cor,
                      SYMBOL = dd$mRNAs)
gene_df <- merge(gene_df,gene,by="SYMBOL")

## geneList 三部曲
## 1.获取基因logFC
geneList <- gene_df$logFC
## 2.命名
names(geneList) = gene_df$ENTREZID
## 3.排序很重要
geneList = sort(geneList, decreasing = TRUE)
head(geneList)
library(org.Hs.eg.db)
#GO富集
GO <- gseGO(
  geneList, #
  ont = "BP",# "BP"、"MF"和"CC"或"ALL"
  OrgDb = org.Hs.eg.db,#人类注释基因
  keyType = "ENTREZID",
  eps = 0,
  pvalueCutoff = 0.5,
  pAdjustMethod = "BH",#p值校正方法
)
#KEGG Pathway富集
gseago <- gseKEGG(geneList, organism = "hsa", pvalueCutoff = 0.05)
gsea<- setReadable(gseago, OrgDb=org.Hs.eg.db,keyType = 'ENTREZID')

ego <- gseGO(geneList     =  geneList,
             OrgDb        = org.Hs.eg.db,
             ont          = "BP", pvalueCutoff = 1,  pAdjustMethod = "BH")
#作图
library(stringr)
library(enrichplot)
dotplot(GO,color = "p.adjust",showCategory=10,split='.sign',font.size = 20,label_format=50)+facet_grid(~.sign, scales = "free")+
  scale_color_gradient(low = "red", high = "blue")

gseaplot2(GO,"GO:0034340",base_size = 12,pvalue_table = F) 
ggsave("SIRT1 dotplot.tiff", device = tiff, width = 12, height = 9,units = "in")
# 绘制 "response to type I IFN" 的GO图

gseaplot2(
  GO, #gseaResult object，即GSEA结果
  "GO:0045088",#富集的ID编号
  title = "GO response to type I interferon", #标题
  color = "green",#GSEA线条颜色
  base_size = 12,#基础字体大小
  rel_heights = c(1.5, 0.5, 1),#副图的相对高度
  subplots = 1:3, #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图，subplots=1#只要第一个图
  pvalue_table = T, #是否添加 pvalue table
  ES_geom = "line" #running enrichment score用先还是用点ES_geom = "dot"
)

head(GO@result$Description, 50)

ggsave("SIRT1 response to type I interferon.tiff", device = tiff, width = 16, height = 12,units = "cm")

#grep 返回匹配元素的字符向量#
grep("interferon|IFN|immune", GO@result$Description, value = TRUE)

# 获取对应的值
GO@result[["p.adjust"]][which(GO@result[["ID"]] == "GO:0034340")]
GO@result[["ID"]][which(GO@result[["Description"]] == "response to type I interferon")]
#可视化网络#
library(aPEAR)

pdf("USP7Network.pdf",height = 10,width = 10)
enrichmentNetwork(GO@result[1:100,])
dev.off()
