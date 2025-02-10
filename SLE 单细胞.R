
#读取文件
seurat <- readRDS("GSE135779_seurat.rds")
library(ggplot2)
library(ggpubr)
library(Matrix)
library(SeuratObject)
library(sp)
library(Seurat)
library(Rcpp)
library(harmony)
library(usethis)
library(devtools)
library(tidyverse)
library(stringr)
library(ggrepel)
#获取目录名字
dir = paste0("E:/SIRT1/BBA/scRNA/21941063/BM/",dir("BM/"))
dir
#标记
names(dir) = c(rep('AA',4 ), rep('Healthy', 4))

#批量读取为一个list
scRNAlist <- list()
for (i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
scRNAlist[[i]] <- CreateSeuratObject(counts, project = names(dir)[i], min.cells=1)
}
#names(scRNAlist) <- name
#merge对象
seurat  <-  merge(scRNAlist[[1]],scRNAlist[2:length(scRNAlist)])

#或者一次性全部读取
counts <- Read10X(data.dir = dir[1])
scRNA1 = CreateSeuratObject(counts, min.cells=1)
dim(seurat)   #查看基因数和细胞总数
table(seurat@meta.data$orig.ident)  #查看每个样本的细胞数
##meta.data添加信息
#proj_name <- data.frame(proj_name=rep("demo2",ncol(scRNA1)))
#rownames(proj_name) <- row.names(scRNA1@meta.data)
#scRNA1 <- AddMetaData(scRNA1, proj_name)
##计算线粒体基因比例
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT[-\\.]")
##绘制小提琴图
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
##设置质控标准
seurat <- subset(seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 15)

#标准化
seurat <- NormalizeData(seurat)

#小提琴图基因表达，添加统计学
my_comparisons <- c('SLE', 'Healthy')#list(c('Healthy', 'Lupus'), c('Lupus', 'non-Les'))
VlnPlot(seurat, features = c("SIRT1"), split.by = "orig.ident")&
  stat_compare_means(method="t.test",hide.ns = F,
                     comparisons = my_comparisons,
                     label="p.signif",
                     bracket.size=0.8, 
                     tip.length=0,
                     size=6)
ggsave("SIRT1.TIFF", device = tiff, width = 8, height = 6, units ="in")

VlnPlot(seurat, features = c("SIRT1"), split.by = "orig.ident")&
  theme_bw()&  
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(color = 'black',face = "bold", size = 12),
        axis.text.y = element_text(color = 'black', face = "bold"),
        axis.title.y = element_text(color = 'black', face = "bold", size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color="black",size = 1.2, linetype="solid"),
        panel.spacing = unit(0.12, "cm"),
        plot.title = element_text(hjust = 0.5, face = "bold.italic"),
        legend.position = 'none')&
  stat_compare_means(method="wilcox.test",hide.ns = F,
                     comparisons = c(my_comparisons),
                     label="p.signif",
                     bracket.size=0.8, 
                     tip.length=0,
                     size=6)&
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))&
  scale_fill_manual(values = c("#FF5744","#208A42", "#FCB31A"))

VlnPlot(seurat, features = c("SIRT1"), split.by = "orig.ident")


#降维，聚类
seurat <- FindVariableFeatures(seurat, nfeatures = 3000)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat, npcs = 50)
seurat <- RunTSNE(seurat, dims = 1:10)
TSNEPlot(seurat)

#去批次
library(harmony)
seurat <- RunHarmony(seurat, group.by.vars = "orig.ident", dims.use = 1:20, max.iter.harmony = 50)
seurat <- RunTSNE(seurat, reduction = "harmony", dims = 1:10)
TSNEPlot(seurat, group.by = "orig.ident")
#分群
seurat <- FindNeighbors(seurat, reduction = "harmony", dims = 1:10) 
seurat <- FindClusters(seurat, resolution = 0.1)
DimPlot(seurat, reduction = "tsne", label = TRUE)
#注释
seurat <- JoinLayers(seurat) 
markers_df <- FindMarkers(object = seurat, ident.1 = 1, min.pct = 0.25)
cl_markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log(1.2))
library(dplyr)
clm<-cl_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(clm,"cl_markers.csv",sep=",",row.names=FALSE)


DimPlot(seurat, reduction = "tsne", label = TRUE)
#singleR自动注释
library(SingleR)
ref_data <- get(load("C:\\R\\library\\celldex\\ref_Human_all.RData")) # 加载参考数据
expr_data <- GetAssayData(seurat, slot = "data")
predictions <- SingleR(test = expr_data, ref = ref_data, labels = ref_data$label.main)

seurat@meta.data$SingleR.labels <-predictions$labels
DimPlot(seurat, group.by = c("seurat_clusters", "SingleR.labels"),reduction = "tsne")

#SLE手动注释
gene.list =c("CD3D",'CD8A','IL7R','CD34',
             'CD14','FCGR3A',
             'MS4A1','CD19','XBP1','MZB1','CD79A',
             'CEACAM8','CD177','S100A8',
             'TPSAB1','TPSB2',
             "GATA1", "HBD",
             "HBA1", "HBA2",
             "CD1C", "CLEC10A", 
             'LILRA4','IRF8','IRF7',
             "APOE", "VCAM1", "DCN")
DotPlot(seurat, features = gene.list)+ggplot2:::coord_flip()
new_ident <- setNames(c('Monocytes',
                        'Granulocyte',
                        'Erythroid progenitors',
                        'Erythroid progenitors',
                        'T/NK cell',
                        'Granulocyte',
                        'Pro/pre‐B cell',
                        'B/Plasma cell',
                        "Monocytes",
                        'Erythroid progenitors',
                        "DC",
                        "EC/FB"
),
                      levels(seurat))
seurat <- RenameIdents(seurat, new_ident)

DimPlot(seurat,  label = TRUE, repel = TRUE)
ggsave("分群注释.pdf", path = "E:\\SIRT1\\BBA\\scRNA\\SLE results\\GSE135779", width= 6 , height= 4 , units='in')

markers <- c("SIRT1","IFNB1")
#添加一列，在两组时就可以就可以显示Average expression 
seurat$cell.type_splitby_group=paste0(seurat@active.ident,"_",seurat$orig.ident)
DotPlot(seurat, 
        features = "SIRT1", 
        group.by = "cell.type_splitby_group"
        )+ RotatedAxis()+ggplot2:::coord_flip()
ggsave("SIRT1 dotplot.tiff",path = "E:\\SIRT1\\BBA\\scRNA\\AA results\\纯PM", width= 10 , height= 4 , units='in')

VlnPlot(seurat, 
        features = "JCHAIN",
        split.by = 'orig.ident',
        #y.max = 1,
        pt.size = 0
        ) 
#筛选数据，解决小提琴图没有数据
a <- subset(seurat,SIRT1 > 0 )
VlnPlot(a, 
        features = "SIRT1",
        split.by = 'orig.ident',
        pt.size = 0
) + stat_compare_means(aes(group = a@meta.data$orig.ident), size = 3)
ggsave("SLE SIRT1.pdf",path = "E:\\SIRT1\\BBA\\scRNA\\SLE results\\GSE135779", width= 10 , height= 6 , units='in')

FeaturePlot(seurat, features = c("CD3D"),
            split.by = 'orig.ident')
library(scran)
FeatureScatter(object = subset(seurat,SIRT1 > 0 ), feature1 = "SIRT1", feature2 = "HALLMARK_INTERFERON_ALPHA_RESPONSE_scaled")
saveRDS(seurat, file="AA-BM_seurat.rds")
saveRDS(cl_markers,file="AA_cl_markers.rds")

#AA注释
new_ident <- setNames(c('T-cell',
                        'Macrophage',
                        'B-cell',
                        'NK cell',
                        'Erythroblast',
                        'Neutrophill',
                        'Fibroblast',
                        'B-cell',
                        'Monocyte',
                        'B-cell',
                        'Endothelial cell',
                        'Erythroblast',
                        'DC',
                        'Erythroblast',
                        'DC',
                        'Megakaryocyte',
                        'Endothelial cell'
),
levels(seurat))
seurat <- RenameIdents(seurat, new_ident)

