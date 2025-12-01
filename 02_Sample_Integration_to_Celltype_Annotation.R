# ==== 加载配置和工作环境 ====
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
gc()
setwd("~/Biostudy/ACE_AI")
getwd()
library(qs)
library(SCP)
library(dplyr)
library(Seurat)
library(ggpubr)
library(cowplot)
library(ggplot2)
library(harmony)
library(patchwork)
library(RColorBrewer)

# ==== 读取数据 ====
ctrl=qread("./CTRL_post_QC.qs")
psmd=qread("./PSMD_post_QC.qs")
ace=qread("./ACE_post_QC.qs")

seurat_list=list(ctrl,psmd,ace)
seurat_obj=Reduce(merge,seurat_list)
print(seurat_obj)

# ==== 降维聚类+去除批次效应 ====
seurat_obj[["RNA"]]=JoinLayers(seurat_obj[["RNA"]])
seurat_obj=NormalizeData(seurat_obj,normalization.method = "LogNormalize",
                         scale.factor = median(seurat_obj@meta.data$nCount_RNA))
seurat_obj=FindVariableFeatures(seurat_obj,nfeatures=3000)
# 我们发现在亚群注释的过程中，如果按照上面的结果直接去做细胞注释
# 会有一些奇奇怪怪的细胞高表达非经典神经内分泌marker
# 我们觉得这样的结果不反映我们的生物学现象，而且这些细胞恰好高表达核糖体基因
# 所以我们希望去除核糖体基因和细胞周期基因对降维聚类的影响
# 计算细胞周期评分
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

# safer matching (recommended)
s.genes.use   <- CaseMatch(search = s.genes,   match = rownames(seurat_obj))
g2m.genes.use <- CaseMatch(search = g2m.genes, match = rownames(seurat_obj))

#seurat_obj[["percent_ribo"]]=PercentageFeatureSet(seurat_obj, pattern = "^Rps|^Rpl")
seurat_obj <- CellCycleScoring(seurat_obj,s.features = s.genes.use,g2m.features = g2m.genes.use)

# 回归掉不感兴趣的变量
seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c("S.Score", "G2M.Score","percent_mt"))

# 使用HVG去跑PCA
seurat_obj=RunPCA(seurat_obj)
ElbowPlot(seurat_obj,ndims = 50)
library(scPCselect)
viz <- visualize_pc_selection(
  seurat_obj,
  max_pcs = 50)
# 根据图表选择合适的PC数
# 图1绿色数字提示为建议选择PC
best_pcs=16

seurat_obj=RunHarmony(seurat_obj,group.by.vars="orig.ident")
seurat_obj=RunUMAP(seurat_obj,dims=1:best_pcs,verbose = T,reduction = "harmony")
seurat_obj=FindNeighbors(seurat_obj,dims = 1:best_pcs,reduction = "harmony")
library(clustree)
library(patchwork)
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  seurat_obj <- FindClusters(seurat_obj, resolution = res)
}
p1 <- clustree(seurat_obj, prefix = 'RNA_snn_res.') + coord_flip()
p2 <- DimPlot(seurat_obj, group.by = 'RNA_snn_res.0.1', label = T)
p1 + p2 + plot_layout(widths = c(3, 1))
#人工挑选RNA_snn_res.0.1这个参数
seurat_obj$seurat_clusters <- seurat_obj$RNA_snn_res.0.1
Idents(seurat_obj) <- "seurat_clusters"

# seurat_obj=FindClusters(seurat_obj,resolution = 0.1)
table(seurat_obj@meta.data$seurat_clusters)
plot1=DimPlot(seurat_obj,reduction = "umap",group.by = "orig.ident",label = T)
plot2=DimPlot(seurat_obj,reduction = "umap",group.by = "seurat_clusters",label = T)+NoLegend()
plot1|plot2

#找差异基因
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5,test.use = "MAST")
# marker基因就两个特点，特异性和高表达
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#10差异不明显，改成20 
top20 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
#marker基因汇总的网站cellmarker

library(GPTCelltype)
library(openai)
top10_frame<-as.data.frame(top20)
Sys.setenv(OPENAI_API_KEY='**********')
#Sys.setenv(OPENAI_API_BASE_URL = "*********")
res <- gptcelltype(top10_frame, tissuename = "Mouse skeletal muscle", model = 'gpt-4-turbo')
# 查看并修改res中某一个值
res

to_mouse_style <- function(x) {
  stringr::str_to_title(tolower(x))  # "CD8A" -> "Cd8a"
}

#T NK 4
DotPlot(seurat_obj,features =c("Cd3d", "Cd3e", "Cd3g", "Cd8a", "Cd4", "Ptprc","Gnly", "Nkg7", "Gzmb"),assay='RNA' )
# B 19
DotPlot(seurat_obj,features =to_mouse_style(c("MS4A1", "CD79A", "IGHM", "CD19", "CD22", "IGKC")),assay='RNA' )
# Neutrophils  20
DotPlot(seurat_obj,features =to_mouse_style(c("ITGAM", "CD14", "S100A8", "S100A9", "CSF3R", "CXCL8", "FCGR3B", "CXCR2", "CMTM2")),assay='RNA' )
# Macrophage/monocyte 8 14
DotPlot(seurat_obj,features =to_mouse_style(c("ITGAM", "CSF1R", "ADGRE1", "ITGB1", "CD14", "CD163", "CD68", "CD86", "MRC1", "LYZ", "MARCO", "LYVE1", "IL1A", "IL1B", "SIGLEC1", "CLEC7A", "CLEC10A")),assay='RNA' )
#GMPs 
DotPlot(seurat_obj,features =to_mouse_style(c("MS4A3", "MPO", "ELANE", "PCLAF", "CTSG", "TPSB2", "PRTN3", "TPSAB1", "CPA3")),assay='RNA' )
# endothelial 7,10,18,21,22
DotPlot(seurat_obj,features =to_mouse_style(c("VWF", "ESAM", "PECAM1", "ENG", "CDH5", "KDR", "TEK", "SLC9A3R2", "CLDN5", "AQP1", "ACKR1", "RAMP3", "PLVAP", "EGFL7", "CD34")),assay='RNA' )
# FAPs 2 3 11 12 13 16 21 
DotPlot(seurat_obj,features =to_mouse_style(c("DCN", "PDGFRA", "LY6E", "COL3A1", "GSN", "LY6A","PIEZO2")),assay='RNA' )
# Adipocyte  17
DotPlot(seurat_obj,features =to_mouse_style(c("ADIPOQ", "PPARG", "CEBPA", "G0S2", "LPL", "PLIN1")),assay='RNA' )
# SMCs 5
DotPlot(seurat_obj,features =to_mouse_style(c("MYH11", "PDGFRB", "MYL9", "ACTA2", "CSPG4", "MYLK")),assay='RNA' )
# Mature skeletal muscle  0 1 9 15 
DotPlot(seurat_obj,features =to_mouse_style(c("MYH1", "ACTA1", "MYL1","MYH7","TNNT3","MYH2","TNNT1")),assay='RNA' )
# MuSCs 6
DotPlot(seurat_obj,features =to_mouse_style(c("PAX7", "MYOD1", "SDC4", "VCAM1", "CXCR4", "MYF5", "DES", "APOC1", "APOE")),assay='RNA' )
# Schwann Cell
DotPlot(seurat_obj,features =to_mouse_style(c("PLP1","GFAP","S100B","MPZ","L1CAM","GAP43","MSLN","ITLN1")),assay='RNA' )
# Pericyte cell
DotPlot(seurat_obj,features =to_mouse_style(c("ACTA2","HBA1","ABCC9","MYH11","DLK","NG2","MYO1B","EGFAM","STEAP4","RGS5","SORBS2")),assay='RNA' )
#Erythrocytes红细胞
DotPlot(seurat_obj,features =c("Hbb-bt","Hbb-bs","Hba1","Cd244","Cd235a","Cd164","Klf1","Cd44","Hbg2","Cd36"),assay='RNA' )
DotPlot(seurat_obj,features =c("Hba1","Hba2","Hbb","Hbd","Hbe1","Hbg1","Hbg2","Hbm","Hbq1","Hbz"),assay='RNA' )
#Plasma cells
DotPlot(seurat_obj,features =c("Mpeg1","Ctss","Ctsb"),assay='RNA' )

# ==== 细胞注释 ====
known_markers <- list(
  "T/NK"       = c("Cd3e", "Cd4", "Cd8a"),
  B            = c("Cd79a", "Ms4a1"),
  Plasma       = c("Mzb1"),
  Epithelium   = c("Epcam", "Krt8", "Krt18"),
  Tumor        = c("Mki67", "Foxj1", "Sox2", "Sox9"),
  "Mono/Macro" = c("Cd14", "Cd68", "Cd163"),
  Mast         = c("Kit"),
  Fibroblast   = c("Col5a2", "Pdgfrb"),
  Pericyte     = c("Cspg4", "Rgs5"),
  Endothelium  = c("Pecam1")
)

known_markers <- list(
  "T/NK" = c(
    "Cd3d", "Cd3e", "Cd3g", "Cd8a", "Cd4", "Ptprc",
    "Gnly", "Nkg7", "Gzmb"
  ),
  
  B = c(
    "Ms4a1", "Cd79a", "Ighm", "Cd19", "Cd22", "Igkc"
  ),
  
  Neutrophil = c(
    "S100a8", "S100a9", "Csf3r",
    "Cxcl8", "Fcgr3b", "Cxcr2", "Cmtm2"
  ),
  
  "Mono/Macro" = c(
    "Itgam", "Csf1r", "Adgre1", "Itgb1", "Cd14",
    "Cd163", "Cd68", "Cd86", "Mrc1", "Lyz",
    "Marco", "Lyve1", "Il1a", "Il1b",
    "Siglec1", "Clec7a", "Clec10a"
  ),
  
  GMP = c(
    "Ms4a3", "Mpo", "Elane", "Pclaf", "Ctsg",
    "Tpsb2", "Prtn3", "Tpsab1", "Cpa3"
  ),
  
  Endothelial = c(
    "Vwf", "Esam", "Pecam1", "Eng", "Cdh5",
    "Kdr", "Tek", "Slc9a3r2", "Cldn5", "Aqp1",
    "Ackr1", "Ramp3", "Plvap", "Egfl7", "Cd34"
  ),
  
  FAP = c(
    "Dcn", "Pdgfra", "Ly6e", "Col3a1", "Gsn",
    "Ly6a", "Piezo2"
  ),
  
  Adipocyte = c(
    "Adipoq", "Pparg", "Cebpa", "G0s2", "Lpl", "Plin1"
  ),
  
  SMC = c(
    "Myh11", "Pdgfrb", "Myl9", "Acta2", "Cspg4", "Mylk"
  ),
  
  SkMature = c(
    "Myh1", "Acta1", "Myl1", "Myh7",
    "Tnnt3", "Myh2", "Tnnt1"
  ),
  
  MuSC = c(
    "Pax7", "Myod1", "Sdc4", "Vcam1", "Cxcr4",
    "Myf5", "Des", "Apoc1", "Apoe"
  ),
  
  Schwann = c(
    "Plp1", "Gfap", "S100b", "Mpz",
    "L1cam", "Gap43", "Msln", "Itln1"
  ),
  
  Pericyte = c(
    "Hba1", "Abcc9", "Dlk",
    "Ng2", "Myo1b", "Egfam", "Steap4", "Rgs5", "Sorbs2"
  )
)

# 经典气泡图，看不同cluster的marker gene表达情况
plot3=DotPlot(object = seurat_obj,
              features = known_markers,
              scale=T,
              group.by = "seurat_clusters")+
  scale_color_gradientn(colors=brewer.pal(9,"Blues"))+
  theme_pubr()+
  theme(axis.text.x = element_text(angle=90)) & NoLegend()


plot3=DotPlot(object = seurat_obj,
              features = known_markers,
              scale=T,
              group.by = "seurat_clusters")+
  scale_color_gradientn(colors=brewer.pal(9,"Blues"))+
  theme_pubr()+
  theme(axis.text.x = element_text(angle=90)) & NoLegend()
plot2|plot3
# 邱同学风格的细胞注释
meta_supp = data.frame(seurat_cluster = 0:(length(unique(seurat_obj$seurat_clusters)) - 1), celltype = NA)
meta_supp[meta_supp$seurat_cluster %in% c(0), 'celltype'] = 'Neutrophils' # "S100a8", "S100a9", "Csf3r"
meta_supp[meta_supp$seurat_cluster %in% c(1), 'celltype'] = 'Mono/Macro' # "Csf1r", "Adgre1","Cd68
meta_supp[meta_supp$seurat_cluster %in% c(2), 'celltype'] = 'B' #"Cd79a", "Ighm", "Cd19"
meta_supp[meta_supp$seurat_cluster %in% c(3), 'celltype'] = 'GMP' # "Ms4a3", "Mpo", "Elane"
meta_supp[meta_supp$seurat_cluster %in% c(4), 'celltype'] = 'Erythroid' #"Hbb-bt","Hbb-bs","Klf1
meta_supp[meta_supp$seurat_cluster %in% c(5), 'celltype'] = 'Endothelial' #"Vwf", "Esam", "Pecam1"
meta_supp[meta_supp$seurat_cluster %in% c(6), 'celltype'] = 'FAP' #"Dcn", "Pdgfra",  "Col3a1"
meta_supp[meta_supp$seurat_cluster %in% c(7), 'celltype'] = 'T/NK' # "Cd3d", "Cd3e", "Cd3g"
meta_supp[meta_supp$seurat_cluster %in% c(8), 'celltype'] = 'B' #"Cd79a", "Ighm", "Cd19"
meta_supp[meta_supp$seurat_cluster %in% c(9), 'celltype'] = 'Mast' #'Tpsb2', 'Ms4a2'
meta_supp[meta_supp$seurat_cluster %in% c(10), 'celltype'] = 'MuSC'  #"Pax7", "Myod1","Myf5
meta_supp[meta_supp$seurat_cluster %in% c(11), 'celltype'] = 'Schwann cells' #"Plp1", "S100b", "Mpz"

for (i in 1:nrow(meta_supp)) {
  seurat_obj@meta.data[which(seurat_obj$seurat_clusters == meta_supp$seurat_cluster[i]), 'celltype_major'] = meta_supp$celltype[i]
}
Idents(seurat_obj) <- 'celltype_major'

# 看看注释情况
plot4=CellDimPlot(
  seurat_obj,
  group.by = "celltype_major",
  theme_use = "theme_blank",
  xlab = "UMAP1",
  ylab = "UMAP2",
  label = TRUE,           
  label_insitu = TRUE,
  show_stat = F,      
  legend.position = "none" 
)
plot5=DotPlot(object = seurat_obj,
              features = known_markers,
              scale=T,
              group.by = "celltype_major")+
  scale_color_gradientn(colors=brewer.pal(9,"Blues"))+
  theme_pubr()+
  theme(axis.text.x = element_text(angle=90)) & NoLegend()
plot4|plot5
table(seurat_obj@meta.data$celltype_major)

DimPlot(seurat_obj , reduction = "umap",label = T) 
DimPlot(seurat_obj, reduction = "umap", split.by ='orig.ident')
DimPlot(seurat_obj, reduction = "umap", group.by='celltype_major')

# ==== 不要忘记保存注释后的对象 ====
qsave(seurat_obj, file = 'seurat_obj_post_annotation.qs')

# ==== 这里顺带把不同种类的细胞也subset了 ====
# 1. 提取所有大类细胞类型
celltypes <- unique(seurat_obj$celltype_major)

# 2. 循环 subset 并分别保存
for (ct in celltypes) {
  # subset 当前这个细胞类型
  obj_sub <- subset(seurat_obj, subset = celltype_major == ct)
  
  # 处理一下名字，避免 /、空格 等在文件名里出问题
  safe_name <- gsub("[^A-Za-z0-9]+", "_", ct)
  
  # 若想在当前目录保存：
  qsave(obj_sub, file = paste0(safe_name, ".qs"))
  
  # 如果你想在一个特定文件夹里保存，例如 "subsets" 目录：
  # dir.create("subsets", showWarnings = FALSE)
  # qsave(obj_sub, file = file.path("subsets", paste0(safe_name, ".qs")))
  
  # 如果你还想把对象放到当前 R 环境里（像 T_NK 那样有变量名），可以加：
  # assign(safe_name, obj_sub)
}
