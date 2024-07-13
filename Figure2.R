library(dplyr)
library(ggplot2)
library(FigR)
library(SummarizedExperiment)
library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10)
library(doParallel)
library(Seurat)
library(Signac)
library(FNN)
my36colors <- c('#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#E5D2DD', '#AB3282', '#BD956A', 
                '#8C549C', '#E0D4CA', '#C5DEBA', '#58A4C3', '#E4C755', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', 
                '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963','#968175')

MCortex <- readRDS('./Mcortex.rds')

#### Figure1A ####
DefaultAssay(Mcortex) <- 'RNA'
p1=DimPlot_scCustom(Mcortex, group.by = "orig.ident", label.size = 4,colors_use = c('#D6E7A3', '#57C3F3'),figure_plot = TRUE)
p2=DimPlot(Mcortex, group.by = "rna_celltype", label = TRUE, repel = TRUE, label.size = 5,cols = my36colors) + 
        theme_dr(xlength = 0.2, ylength = 0.2, arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
        theme(panel.grid = element_blank(), axis.title = element_text(face = 2,hjust = 0.03)) + ggtitle("Gene expression")
DefaultAssay(Mcortex) <- 'ATAC'
p3=DimPlot_scCustom(Mcortex, group.by = "orig.ident", label.size = 4,colors_use = c('#D6E7A3', '#57C3F3'),figure_plot = TRUE)
p4=DimPlot(Mcortex, group.by = "rna_celltype", label = TRUE, repel = TRUE, label.size = 5,cols = my36colors) + 
        theme_dr(xlength = 0.2, ylength = 0.2, arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
        theme(panel.grid = element_blank(), axis.title = element_text(face = 2,hjust = 0.03)) +  ggtitle("Chromatin accessibilty")
(p1|p2)/(p3|p4)

#### Figure1B ####
DefaultAssay(Mcortex) <- 'RNA'
options(repr.plot.width=18, repr.plot.height=10)
makerlist1<- c('Hap1','Myo16','Nectin3','Tle4','Tac1','Ndst4','Itga8','Lhx6','Pvalb','Inpp4b','Slc1a2','Slc1a3','Lama2',
               'Mobp','Mbp','Mog','Mag','Plp1','Satb2','Reln','Cnr1','Sox6','Pdgfra','Brinp3','Tnr','Fras1',
               'Shisa6','Apbb1ip','P2ry12','Dock8','Hexb','Flt1','Vwf','Pecam1','Rgs6','Th','Slc6a3')
Idents(object = Mcortex) <- "rna_celltype"
p5 <- Seurat::DotPlot(Mcortex, features = makerlist1, group.by='rna_celltype', scale = TRUE, cluster.idents=FALSE, col.min = 0,col.max = 3, dot.scale=9) + RotatedAxis() + 
     scale_color_gradientn(colours = c(rgb(252/255, 247/255, 243/255),rgb(253/255, 200/255, 180/255),rgb(229/255, 51/255, 38/255),rgb(102/255, 0/255, 13/255))) + #(colours = c(rgb(252/255, 247/255, 243/255),rgb(237/255, 213/255, 192/255),rgb(224/255, 147/255, 121/255),rgb(176/255, 49/255, 51/255),rgb(81/255, 24/255, 29/255)))
     xlab('marker genes') + ylab('Celltypes') + scale_y_discrete(limits=rev(levels(Mcortex))) + scale_x_discrete(limits=makerlist1)
p5

#### Figure1C ####
DefaultAssay(Mcortex) <- 'ATAC'
genes.use=c('Nefh','Slc17a7','Gad1','Vip','Slc1a3','Plp1','Pdgfra')
# link peaks to genes
Mcortex <- LinkPeaks(object = Mcortex, peak.assay = "ATAC", expression.assay = "RNA", genes.use = genes.use)
Idents(Mcortex) <- 'rna_celltype'
DefaultAssay(Mcortex) <- 'ATAC'
options(repr.plot.width=25, repr.plot.height=35)
cov_plot <- CoveragePlot(object = Mcortex,region = genes.use,annotation = TRUE, peaks = FALSE, links = FALSE)  & scale_fill_manual(values = my36colors)
cov_plot

#### Figure1D ####
# DORCs score matrix was processed according to FigR documentation
# SCARlink score matrix was processed according to SCARlink documentation
# gene activity was calculated using Signac
MCortex_10x <- readRDS('./Mcortex_10x.rds')
MCortex_scCAT <- readRDS('./Mcortex_scCAT.rds')
DefaultAssay(cortex_scCAT) <- 'RNA'
DefaultAssay(cortex_10x) <- 'RNA'
cortex_scCAT=FindVariableFeatures(cortex_scCAT, nfeatures = 3000)
cortex_10x=FindVariableFeatures(cortex_10x, nfeatures = 3000)
hvg_scCAT=VariableFeatures(cortex_scCAT[['RNA']])
hvg_10x=VariableFeatures(cortex_10x[['RNA']])
#
rnamat_scCAT <- cortex_scCAT@assays$RNA$data
rnamat_10x <- cortex_10x@assays$RNA$dat
dorcmat_10x <- readRDS('./DORCmat_10x.rds') # or SCARlinkmat geneactivity
dorcmat_scCAT <- readRDS('./DORCmat_scCAT.rds') # or SCARlinkmat geneactivity
# C4
overlap_gene_scCAT=intersect(hvg_scCAT,rownames(dorcmat_scCAT))
rnamat_scCAT=rnamat_scCAT[overlap_gene_scCAT,]
dorcmat_scCAT=dorcmat_scCAT[overlap_gene_scCAT,]
# 10X
overlap_gene_10x=intersect(hvg_10x,rownames(dorcmat_10x))
rnamat_10x=rnamat_10x[overlap_gene_10x,]
dorcmat_10x=dorcmat_10x[overlap_gene_10x,]
# corr compute 
corr_compute <- function(rnamat, dorcmat, overlapgenes, used_method='') {
    # Initialize empty lists/vectors to store results
    corrs <- c()
    gene_name <- c()
    # Loop over each gene in overlapgenes
    for (gene in overlapgenes) {
        # Extract RNA and ATAC matrices based on the gene and calculate correlation
        rna <- as.numeric(rnamat[gene,])
        atac <- as.numeric(dorcmat[gene,])
        m_obs <- cor.test(rna, atac, method = "spearman")$estimate
        corrs <- c(corrs, m_obs)
        # Append gene name to list
        gene_name <- c(gene_name, gene)
    }
    # Create a data frame to store correlations and gene names
    df <- data.frame(corr = corrs, gene = gene_name, method = used_method)
    return(df)
}
corr_scCAT <- corr_compute(rnamat_scCAT, dorcmat_scCAT, rownames(dorcmat_scCAT), used_method='DORC_C4')
corr_10x <- corr_compute(rnamat_10x, dorcmat_10x, rownames(dorcmat_10x), used_method='DORC_10x')
# corr plot
options(repr.plot.width=6, repr.plot.height=5.4, repr.plot.res = 150)
p1=ggplot(data = SCARlink, aes(x = method, y = corr, fill = method)) + geom_boxplot(width = 0.6, alpha = 0.5,outlier.size = 0.4) + ylim(-0.3,0.8) +
        scale_fill_manual(values=c("DORC_10x"="black","DORC_scCAT"="red",'Geneactivity_10x'='black','Geneactivity_scCAT'='red',"SCARlink_10x" = "black", "SCARlink_scCAT" = "red")) + theme(legend.position='none') +
        theme(panel.background = element_rect(fill = "transparent",colour = NA),panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),axis.text.x = element_text(angle = 45, hjust = 1))
#geom_errorbar(aes(ymin = corr - error, ymax = corr + error), width = 0.2, position = position_dodge(0.9))
p1

#### Figure1E ####
umap.C4 <- as.data.frame(C4@reductions[["umap"]]@cell.embeddings)
umap.10x <- as.data.frame(mul@reductions[["umap"]]@cell.embeddings)
# add smoothed matrix
gene='Tle4'
options(repr.plot.width=20, repr.plot.height=10, repr.plot.res = 150)
p1 <- plotMarker2D(umap.10x,RNAmat.s_10x,markers = c(gene),pointSize = 0.1,maxCutoff = "q0.99",colorPalette = "brewer_purple") + ggtitle("10X_RNA")
p2 <- plotMarker2D(umap.C4,RNAmat.s_C4,markers = c(gene),pointSize = 0.1,maxCutoff = "q0.99",colorPalette = "brewer_purple") + ggtitle("scCAT_RNA")
p3 <- plotMarker2D(umap.10x,geneactivity_mat_10x,markers = c(gene), pointSize = 0.1,maxCutoff = "q0.99",colorPalette = "brewer_heat") + ggtitle("10X_gene_activity")
p4 <- plotMarker2D(umap.C4,geneactivity_mat_C4,markers = c(gene),pointSize = 0.1,maxCutoff = "q0.99",colorPalette = "brewer_heat") + ggtitle("scCAT_gene_activity")
p5 <- plotMarker2D(umap.10x,dorcMat.s_10x,markers = c(gene), pointSize = 0.1,maxCutoff = "q0.99",colorPalette = "brewer_heat") + ggtitle("10X_DORCs")
p6 <- plotMarker2D(umap.C4,dorcMat.s_C4,markers = c(gene),pointSize = 0.1,maxCutoff = "q0.99",colorPalette = "brewer_heat") + ggtitle("scCAT_DORCs")
p7 <- plotMarker2D(umap.10x_d,SCARlink_10x.s,markers = c(gene),pointSize = 0.1,maxCutoff = "q0.99",colorPalette = "brewer_heat") + ggtitle("scCAT_SCARlink")
p8 <- plotMarker2D(umap.C4_d,SCARlink_C4.s,markers = c(gene),pointSize = 0.1,maxCutoff = "q0.99",colorPalette = "brewer_heat") + ggtitle("scCAT_SCARlink")
# (p1|p3|p5)/(p2|p4|p6)
(p1|p3|p5|p7)/(p2|p4|p6|p8)