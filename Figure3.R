library(dplyr)
library(ggplot2)
library(Seurat)
library(TFBSTools)
library(presto)
library(motifmatchr)
library(qlcMatrix)
library(scCustomize)
library(RColorBrewer)
library(ggrepel)
library(GenomicRanges)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(chromVAR)
library(pheatmap)
library(zoo)
library(ComplexHeatmap)
library(cowplot)

#### Figure3A ####
Mskin <- readRDS('./Mskin.rds')
Mskin <- FindMultiModalNeighbors(seu, reduction.list = list("harmony", "integrated_lsi"), dims.list = list(1:40, 2:40))
Mskin <- RunUMAP(Mskin, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
Mskin <- FindClusters(Mskin, graph.name = "wsnn", algorithm = 3, resolution = c(0.1, 0.3, 0.5), verbose = FALSE)
cluster_cols = c('#6778AE','#F7F398','#91D0BE','#B2DF8A','#58A4C3','#aa8282','#666666','#476D87','#55A1B1','#E39A35','#DCC1DD','#E5D2DD','#A6761D','#CCE0F5','#BEAED4')
p1 <- DimPlot(Mskin, group.by = "rna_celltype", label = TRUE, repel = TRUE,pt.size=0.3, label.size = 5,cols = cluster_cols, reduction = "wnn.umap") + 
        theme_dr(xlength = 0.2, ylength = 0.2, arrow = arrow(length = unit(0.2, "inches"),type = "closed")) + 
        theme(panel.grid = element_blank(), axis.title = element_text(face = 2,hjust = 0.03)) + ggtitle("WNN UMAP")
DefaultAssay(Mcortex) <- 'RNA'
p2 <- DimPlot(Mcortex, group.by = "rna_celltype", label = TRUE, repel = TRUE, label.size = 5,cols = cluster_cols, reduction = "umap") + 
        theme_dr(xlength = 0.2, ylength = 0.2, arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
        theme(panel.grid = element_blank(), axis.title = element_text(face = 2,hjust = 0.03)) + ggtitle("Gene expression")
DefaultAssay(Mcortex) <- 'ATAC'
p3 <- DimPlot(Mcortex, group.by = "rna_celltype", label = TRUE, repel = TRUE, label.size = 5,cols = cluster_cols,  reduction = "atac.umap") + 
        theme_dr(xlength = 0.2, ylength = 0.2, arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
        theme(panel.grid = element_blank(), axis.title = element_text(face = 2,hjust = 0.03)) +  ggtitle("Chromatin accessibilty")
p1|p2|p3

#### Figure3D ####
# Find DEGs bw clusters
DefaultAssay(seu) <- "SCT"
Idents(seu) <- seu$rna_celltype
de_genes <- FindAllMarkers(seu,only.pos = T,logfc.threshold = 0.1,max.cells.per.ident = 1000)
deg <- unique(de_genes[de_genes$p_val_adj<0.05 & de_genes$avg_log2FC>0.1,'gene'])
# gene-peak linkage
DefaultAssay(seu) <- "ATAC"
link <- Links(seu)
link <- data.frame(link)
link$adj_pval <- p.adjust(link$pvalue,method = 'BH')
link$gene_cluster <- de_genes[match(link$gene,de_genes$gene),'cluster'] 
lvl <- levels(seu)
link <- arrange(link,factor(gene_cluster,levels=lvl))
link <- link[link$adj_pval<0.05 & link$score>0,]
t <- link %>% group_by(gene) %>% summarise(n=n())
summary(t$n) 

# heatmap (gene expression)
link <- link[link$score>0.1,]
link <- arrange(link,factor(gene_cluster,levels=lvl))
DefaultAssay(seu) <- "RNA"
gene.use <- unique(link$gene)
cc.genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
seu <- CellCycleScoring(seu, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE) 
seu <- ScaleData(seu,vars.to.regress = c("S.Score", "G2M.Score"),features=gene.use)
rna <- seu@assays$RNA@scale.data[gene.use,]
asplit_cells <- split(colnames(seu), seu$rna_celltype)
n <- 10
means <- do.call(cbind, lapply(lvl, function(x){
        df <- rna[,asplit_cells[[x]]]
        t(apply(df,1,function(x){rollapply(x,n,mean,by=n)}))
}))
celltype_major <- unlist(lapply(lvl, function(x) rep(x, length(asplit_cells[[x]])%/%n)))
anno_col <- data.frame(celltype_major)
rownames(anno_col) <- colnames(means) <- paste0(seq(1:ncol(means)),colnames(means))
anno_col <- data.frame(celltype_major)
rownames(anno_col) <- colnames(means) <- paste0(seq(1:ncol(means)),colnames(means))
# 
lvl <- c('Basal','Spinous','Fibroblasts','Chod','Dermal_papilla','Cycling','Vascular_endothelium','Lymphatic_endothelium','Muscle','Melanocytes','Schwann','Mast_cells','Mac','Neural_crest','Periderm')
anno_col <- arrange(anno_col,factor(celltype_major,levels=lvl))
anno_col$celltype_major <- factor(anno_col$celltype_major,levels=lvl)
# mark gene
genes_t <- c("Krt5",'Col17a1','Krt15','Krt1','Krt10','Col1a1','Lum','Col1a2','Col2a1','Col9a3','Acan','Cdh5','Flt4',
           'Lyve1','Prox1','Rgs5','Mitf','Sox10','Pax3','Itgb8','Scn7a','Kit','Tpsb2','Cpa3','Il1rl1','C1qb','Krt17','Grhl3')
genes <- as.data.frame(genes_t)
first_occurrence <- sapply(genes$genes, function(x) which(rownames(test) == x)[1])
BluYl = colorRampPalette(colors = c("blue","yellow"))(100)
label_colors <- list(celltype_major=c('Fibroblasts' = '#91D0BE','Dermal_papilla' = '#58A4C3','Periderm' = '#B2DF8A','Basal' = '#6778AE',
  'Spinous' = '#F7F398','Neural_crest' = '#CCE0F5','Melanocytes' = '#E39A35','Schwann' = '#DCC1DD','Mast_cells' = '#E5D2DD','Mac' = '#A6761D',
  'Vascular_endothelium' = '#666666','Lymphatic_endothelium' = '#476D87','Muscle' = '#55A1B1','Cycling' = '#aa8282','Chod' = '#625D9E'))
p1 = pheatmap(means[link$gene,rownames(anno_col)], cluster_rows = F, cluster_cols = F, scale = "row", show_rownames = F,
         breaks = seq(-2, 2, length = 101), col = BluYl,
         annotation_col = anno_col,show_colnames = F, annotation_colors = label_colors, use_raster = TRUE,
         main='Gene expression') + rowAnnotation(link = anno_mark(at = first_occurrence, labels = genes$genes_t, labels_gp = gpar(fontsize = 8)))

# heatmap (chromatin accessibility)
peak.use <- unique(link$peak)
DefaultAssay(seu) <- "ATAC"
Idents(seu) <- 'rna_celltype'
atac <- RunTFIDF(seu@assays$ATAC@counts[peak.use,])
asplit_cells <- split(colnames(seu), seu$rna_celltype)
means_peak <- do.call(cbind, lapply(lvl, function(x){
        df <- atac[,asplit_cells[[x]]]
        t(apply(df,1,function(x){rollapply(x,n,mean,by=n)}))
}))
celltype_major <- unlist(lapply(lvl, function(x) rep(x, length(asplit_cells[[x]])%/%n)))
anno_col2 <- data.frame(celltype_major)
rownames(anno_col2) <- colnames(means_peak) <- paste0(seq(1:ncol(means_peak)),colnames(means_peak))
anno_col2$celltype_major <- factor(anno_col2$celltype_major,levels=lvl)
anno_col2 <- arrange(anno_col2,factor(celltype_major,levels=lvl))
options(repr.plot.width = 12, repr.plot.height = 9)
p2 = pheatmap(means_peak[link$peak,rownames(anno_col2)],cluster_rows = F, cluster_cols = F, scale = "row", show_rownames = F,
               breaks = seq(-2, 2, length = 101), color = hcl.colors(100, "BluYl"),
               annotation_col = anno_col, show_colnames = F, annotation_colors = label_colors, use_raster=TRUE,
               main='Chromatin accessibility') + rowAnnotation(link = anno_mark(at = first_occurrence, labels = genes$genes_t, labels_gp = gpar(fontsize = 8)))
                           
                           
#### Figure3D ####
DefaultAssay(Mskin) <- "ATAC"
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(x = JASPAR2020,opts = list(species = 9606, all_versions = FALSE))
# add motif information
Mskin <- AddMotifs(object = Mskin,genome = BSgenome.Mmusculus.UCSC.mm10,pfm = pfm)
Mskin <- RunChromVAR(object = Mskin, genome = BSgenome.Mmusculus.UCSC.mm10)
DefaultAssay(seu) <- 'RNA'
p1 <- FeaturePlot_scCustom(seurat_object = seu, colors_use = viridis_light_high, features = "Twist2", min.cutoff=0, max.cutoff=3)
DefaultAssay(seu) <- 'aRNA'
p2 <- FeaturePlot_scCustom(seurat_object = seu, colors_use = viridis_magma_dark_high, reduction='umap.atac', features = "Twist2", min.cutoff=0, max.cutoff=1.5)
DefaultAssay(seu) <- 'chromvar'
p3 <- FeaturePlot_scCustom(seurat_object = seu, colors_use = c("grey","red"), reduction='umap.atac', features = "MA0633.1", min.cutoff=0, max.cutoff=6)
DefaultAssay(seu) <- 'ATAC'
p4 <- MotifPlot(object = seu, motifs = "MA0633.1") 
plot_grid(p1, p2, p3, p4, axis = "b", align = "h",ncol = 4)