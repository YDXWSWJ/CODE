# ==============================================================

#!!!!!!R 4.2.3  跟R4.0.2画出的图还不太一样！！！！！！！！！！
# 据判断 R4.2.3应该是准的
# ISG score
# 先测试下IFN-a pathway score
# AUcell
library(AUCell)
library(scales)
library(clusterProfiler)
# create list

# AUCell
countexp <- seurat@assays$RNA@counts
# step1:排序
cells_rankings <- AUCell_buildRankings(countexp,nCores = 1, 
                                       plotStats = F)
gmtFile <- "E:/SIRT1/BBA/scRNA/HALLMARK_INTERFERON_ALPHA_RESPONSE.v2024.1.Hs.gmt"
geneSets <- GSEABase::getGmt(gmtFile)
#gene<-as.list(IL2_STAT5_pathway)

# step2:计算
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
#aucMaxRank=10)

# step3:结果
signature_exp <- data.frame(getAUC(cells_AUC))
#GYM@assays$METABOLISM$score <- signature_exp
#metabolism.matrix <- GYM@assays$METABOLISM$score
#metabolism.matrix <- signature_exp
matrix <- t(signature_exp)
metadata <- seurat@meta.data
metadata <- cbind(metadata,matrix)
seurat@meta.data <- metadata
#summary(scRNA@meta.data$gene)
#Fcols = c('#FFFF00','#FFD700','#FF8C00','#FF4500','#FF0000','#A52A2A','#8B0000','#800000','#000000')
#FeaturePlot(scRNA,features = c('gene'),cols =Fcols)
colnames(seurat@meta.data)

seurat$HALLMARK_INTERFERON_ALPHA_RESPONSE_scaled <- rescale(seurat$HALLMARK_INTERFERON_ALPHA_RESPONSE,
                                                           to = c(0,10))

Fcol <- c('#FFFF00','#FFD700','#FF8C00','#FF4500','#FF0000','#A52A2A','#8B0000','#800000','#000000')


Fcol1 <- c('#00008B','#4169E1','#006400','#90EE90','#FFD700','#FFA500','#8B0000','#800000','#000000')
Fcol1 <- c('#00008B','#4169E1','#FFD700','#FFA500','#8B0000','#800000','#000000')


Fcol2 <- c('#494c85','#6fbf6b','#e5e138','#ff8c0a')


# windowsFonts(HEL = windowsFont("Helvetica CE 55 Roman"),
#              RMN = windowsFont("Times New Roman"),
#              ARL = windowsFont("Arial"))

seurat$score_of_interferon_alpha_response <- seurat$HALLMARK_INTERFERON_ALPHA_RESPONSE_scaled
library(ggplot2)
p <- FeaturePlot(seurat,
                 features = c("HALLMARK_INTERFERON_ALPHA_RESPONSE_scaled"),
                 split.by = 'orig.ident',
                 cols = Fcol1,
                 pt.size = 0.5)+
  #ggtitle('score') +
  theme(legend.text = element_text(size = 10),
        legend.position = 'right')

p

saveRDS(RMacro_c,file = 'RMacro_scRNA_ex_068_harmony_c_206.rds')
pdf(file = "score_IFNa_response_lableT.pdf", width = 12,height = 10)
p
dev.off()


png(filename = "featureplot_IFNa_score2.png", width = 8000,height = 5000, units = "px", res = 600)
p
dev.off()


VlnPlot(seurat,features = c("HALLMARK_INTERFERON_ALPHA_RESPONSE_scaled"),
             split.by = 'orig.ident',
             pt.size = 0) + stat_compare_means(aes(group = seurat@meta.data$orig.ident), size = 3)
p
ggsave(p,filename = "VlnPlot_IFNa_score.png",width =15,height = 8,path = savepath)
saveRDS(scRNA,file = 'scRNA_206_harmony.rds')
