# Script used to generate the figures - ccRCC scRNAseq dataset example
# 8/11/2023
# Lucile Massenet-Regad

# Load libraries
library(BiocGenerics)
library(ggplot2)
library(dplyr)
library(icellnet)
library(gridExtra)
library(Seurat)
library(ComplexHeatmap)
library(circlize)


# Load ICELLNET database
db=as.data.frame(readxl::read_excel("ICELLNET/Databases/ICELLNETdb_V2.xlsx", sheet = 1, guess_max = 3000))
db.name.couple=name.lr.couple(db)
assay="RNA" #SCT or RNA

# load seurat object, retrieve matrix, average expression for each cluster
seurat=readRDS("~/Downloads/GSE222703_SeuratObj_ALL_integrated.rds")

seurat$Celltype_Harmony2[which(seurat$Celltype_Harmony2=="Mac_C1QA")]="Macrop" #rename macrophage cluster
seurat$Celltype_Harmony2[which(seurat$Celltype_Harmony2=="ccRCC2")]="Tumor_cell" #rename ccRCC2 as tumor cell cluster
Idents(seurat)=seurat$Celltype_Harmony2

seurat.tum = subset(seurat, cells=which(seurat$Tissue=="Tum"))
rm(seurat)

Idents(seurat.tum)=seurat.tum$Celltype_Harmony2
table(Idents(seurat.tum))

data <- sc.data.cleaning(object = seurat.tum, db = db, filter.perc = 10, save_file = T, force.file = F, path= paste0(work.dir,output.dir,'ICELLNET_', assay, "/"))

# Gene scaling for ICELLNET
data.icell= as.data.frame(gene.scaling(data, n=1, db=db))

# central cell and partner cell selection
target=data.frame("ID"=colnames(data.icell), "Cell_type"=colnames(data.icell), "Class"=colnames(data.icell))
rownames(target)=target$ID
CC="Tumor_cell"
PC=c("Treg","CD8T", "CD4T", "NK", "Macrop","cDC2", "Endoth","Fibro")


#direction communication score
direction="out"

# data selection and ICELLNET score
CC.data=as.data.frame(data.icell[,c(CC, "Symbol")], row.names = rownames(data.icell))
PC.data=as.data.frame(data.icell[,c(PC, "Symbol")], row.names = rownames(data.icell))


score.computation.2= icellnet.score(direction=direction, PC.data=PC.data, 
                                    CC.data= CC.data,  
                                    PC.target = target, PC=PC, CC.type = "RNAseq", 
                                    PC.type = "RNAseq",  db = db)
score2=as.data.frame(score.computation.2[[1]])
lr2=score.computation.2[[2]]
ymax=max(score2)


# Figure C
PC.col=c("B_cell"="lightblue", "Treg"="lightblue","CD8T"="lightblue", "CD4T"="lightblue", "NK"="lightblue", "Macrop"="white","cDC2"="white", "Neutrop"="grey", "Endoth"="salmon","Fibro"="salmon", "Tumor_cell"="salmon")
network.create(icn.score = score2, scale = c(round(min(score2-1)),round(max(score2+1))), direction = "out", PC.col=PC.col)

# Figure D
LR.family.score(lr=lr2,db.couple=db.name.couple, plot=T, title=paste0(CC,"-", direction),  ymax = ymax, cluster_rows = T, cluster_columns = T) 

# Figure E
LR.heatmap(lr=lr2, thresh = 30, topn = 10, sort.by = "var", db.couple = db.name.couple, title = "10 Most different interactions") 

# Figure F
LR.viz(data = data.icell[,c(PC, CC),], int = "CD70 / CD27", family=NULL, db = db, plot=T)
LR.viz(data = data.icell[,c(PC, CC),], int = "PROS1 / AXL", family=NULL, db = db, plot=T)
