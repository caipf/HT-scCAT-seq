library(Signac)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(viridis)
library(zoo)
library(pheatmap)
library(reshape2)

# load data
skin_fib <- readRDS('./skin_fibroblast.rds')
#### FigureC-G ####
# Load diffusion pseudotime, terminal prob, and absorption rate data
skin_fib_meta <- read.csv('./Mskin_fib_meta.csv',sep='\t',header = TRUE,row.names=1)
colnames(skin_fib_meta) <- c('diffusion_pseudotime','terminal_prob', 'transition_Fib.DC','transition_Fib.Deep','transition_Fib.Inter','transition_Fib.Lower') 
skin_fib@meta.data <- cbind(skin_fib@meta.data, skin_fib_meta) 
# load driver gene
driver_gene <- read.csv('./Mskin_fib_drivergene.csv',header = TRUE)
rownames(driver_gene) = driver_gene$X
driver_DC <- driver_gene[, grep("Fib.DC",colnames(driver_gene))]
driver_DC <- driver_DC[order(driver_DC$Fib.DC_pval), ]
## fib.DC ##
seu$tmp <- seu$transition_Fib.DC > summary(seu$transition_Fib.DC)['Mean'] 
tmp <- seu[,seu$tmp]
corr.thres <- 0.05
Fib.DC <- na.omit(driver_DC[order(driver_DC$Fib.DC_corr,decreasing = T),]) 
Fib.DC$adjusted_pval <- p.adjust(Fib.DC$Fib.DC_pval,method="BH") 
Fib.DC <- Fib.DC[Fib.DC$adjusted_pval < 0.05 & Fib.DC$Fib.DC_corr > corr.thres,]

# heatmap
seu$tmp <- seu$transition_Fib.Deep > summary(seu$transition_Fib.Deep)['Mean']
tmp <- seu[,seu$tmp]
corr.thres <- 0.05
Fib.Deep <- na.omit(driver_Deep[order(driver_Deep$Fib.Deep_corr,decreasing = T),])
Fib.Deep$adjusted_pval <- p.adjust(Fib.Deep$Fib.Deep_pval,method="BH") 
Fib.Deep <- Fib.Deep[Fib.Deep$adjusted_pval < 0.05 & Fib.Deep$Fib.Deep_corr > corr.thres,]    
gene.use <- rownames(Fib.Deep) 
DefaultAssay(tmp) <- "RNA"
tmp <- ScaleData(tmp,vars.to.regress = c("S.Score", "G2M.Score"),features=gene.use) 
pt <- tmp$diffusion_pseudotime
rna <- tmp@assays$RNA@scale.data[gene.use,names(pt)[order(pt)]] 
n <- 100
means <-t(apply(rna,1,function(x){rollapply(x,n,mean,by=n)}))
anno_col <- data.frame(diffusion_pseudotime = rollapply(pt[order(pt)],n,mean,by=n),
                       fate_probability = rollapply(tmp$transition_Fib.Deep[order(pt)],n,mean,by=n),
                       terminal_likelihood = rollapply(tmp$terminal_prob[order(pt)],n,mean,by=n),
                       developmental_stage =  rollapply(tmp$stage[order(pt)],n,mean,by=n))
rownames(anno_col) <- colnames(means) <- paste0(seq(1:ncol(means)),colnames(means))
peaks <- apply(means, 1, function(x) which(x == max(x))) 
rowOrd <- order(peaks)
pheatmap(means[rowOrd,],cluster_rows = F, cluster_cols = F, scale = "row",
         breaks=seq(-2,2,length=101),
         # gaps_row = c(cut1,cut2),
         annotation_col = anno_col,show_colnames = F, show_rownames = F)
               
# cluster               
q <- quantile(anno_col$diffusion_pseudotime,probs = seq(0,1,0.33))
cut1 <- max(which(anno_col$diffusion_pseudotime<=q[2]))
cut1 <- length(peaks[peaks<=cut1])
cut2 <- max(which(anno_col$diffusion_pseudotime<=q[3]))
cut2 <- length(peaks[peaks<=cut2])
genes <- gene.use[rowOrd]
t1 <-genes[1:cut1]
t2 <- genes[cut1:cut2]
t3 <- genes[cut2:length(genes)]
plot_multipleGene <- function(gene.use){
        df <- data.frame(expr=t(means[gene.use,]),diffusion_pseudotime=anno_col$diffusion_pseudotime)
        melted <- melt(df,id.vars = 'diffusion_pseudotime')
        ggplot(melted, aes(y = value, x = diffusion_pseudotime,group=variable)) +
                geom_point()+theme_bw()+
                geom_smooth(method = "loess", span = 1,color='darkgrey') +
                ylab("Relative expression") + xlab("Diffusion pseudotime")
}

p1 <- plot_multipleGene(t1[1:10])+ggtitle("Start")
p2 <- plot_multipleGene(t2[1:10])+ggtitle("Middle")
p3 <- plot_multipleGene(tail(t3,10))+ggtitle("End")
p1+p2+p3

# Enrichment analysis
library("WebGestaltR")
runEnrich <- function(gene,thr,name){
        WebGestaltR(enrichMethod="ORA", organism="mmusculus",
                    enrichDatabase=c("geneontology_Biological_Process_noRedundant"),
                    enrichDatabaseType="genesymbol",
                    interestGene=gene,
                    minNum=5,maxNum = 500,
                    interestGeneType="genesymbol",fdrThr = thr,
                    referenceSet ="genome_protein-coding",
                    referenceGeneType="genesymbol",isOutput = T,
                    #outputDirectory = "GSEA/",
                    projectName=name)  
}
lst <- list(t1,t2,t3)
names(lst) <- c("start","middle","end")
df <- do.call(rbind,lapply(names(lst),function(x){
        write.table(lst[[x]],paste0("GSEA/anterior_",x,".txt"),row.names = F,col.names = F,quote=F,sep="\t")
        a <- runEnrich(gene=lst[[x]],thr=0.05,name=paste0("anterior_",x))
        #if (is.null(a)) {a <- runEnrich(gene=lst[[x]],thr=0.,name=paste0("anterior_",x))}
        a <- a[order(a$enrichmentRatio,decreasing = T),]
        df <- data.frame(a[1:15,c(2,7,9)])
        df$cluster <- x
        df
}
))
df <- na.omit(data.frame(df))
df$cluster <- factor(df$cluster,levels=c("start","middle","end"))
ggplot(df,aes(x=cluster,y = factor(description,levels=rev(unique(description))), 
              color = -log10(FDR), size = enrichmentRatio)) + 
        geom_point() + scale_color_viridis_c(name = '-log10(FDR)') + 
        cowplot::theme_cowplot()+ylab("")+xlab("")
               
#### FigureH-I ####
library(FigR)
library(ComplexHeatmap)
library(scCustomize)
library(JASPAR2020)
library(Signac)
library(Seurat)
library(TFBSTools)
# library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)
library(cowplot)
set.seed(123)
## DORC analysis
ATAC.se = SummarizedExperiment(assays = list(counts=skin_fib@assays$ATAC@counts),rowRanges = skin_fib@assays$ATAC@ranges,colData = skin_fib@meta.data)
RNAmat = skin@assays$RNA@data
# running codes
cisCorr <- FigR::runGenePeakcorr(ATAC.se = ATAC.se,
                           RNAmat = RNAmat,
                           genome = "mm10", # One of hg19, mm10 or hg38 
                           nCores = 8,
                           p.cut = NULL, # Set this to NULL and we can filter later
                           n_bg = 100)
cisCorr.filt <- cisCorr %>% filter(pvalZ <= 0.05)
dorcGenes <- dorcJPlot(dorcTab = cisCorr.filt, cutoff = 5, labelTop = 10, cleanLabels=TRUE, returnGeneList = TRUE, force=20)
# DORC calling
dorcMat <- getDORCScores(ATAC.se = ATAC.se, dorcTab = cisCorr.filt, geneList = dorcGenes, nCores = 4)
# DORC-TF
figR.d <- runFigRGRN(ATAC.se = ATAC.se, dorcTab = cisCorr.filt, genome = "mm10", dorcMat = dorcMat, rnaMat = RNAmat, nCores = 5)
# filter score > 2
res.h <- plotfigRHeatmap(figR.d = figR.d,score.cut = 2,column_names_gp = gpar(fontsize=6),show_row_dend = TRUE)
draw(res.h,heatmap_legend_side = "right",gap = unit(0, "cm"))
require(ggplotify)
g = as.ggplot(res.h)
res.h.mat <- res.h@matrix
set.seed(123);row_idx <- row_order(res.h)
set.seed(123);col_idx <- column_order(res.h)
rowid_dorc_list <- rownames(res.h.mat)[row_idx]
colid_tf_list <- colnames(res.h.mat)[col_idx]
regMat <- res.h.mat[row_idx,col_idx]

# module GRNplot
library(ggraph)
library(tidygraph)
library(tidyverse)
library(igraph)
sliceModule <- function(mat = regMat, row_sel = c('LPL','MPEG1'), col_sel = c('STAT5A','SHOX')) {
    rowid <- rownames(mat)
    colid <- colnames(mat)
    
    rowid_start <- grep(paste0("^",row_sel[1],"$"), rowid)
    rowid_end <- grep(paste0("^",row_sel[2],"$"), rowid)
    stopifnot(rowid_start < rowid_end)
    row_slice <- rowid_start:rowid_end
    
    colid_start <- grep(paste0("^",col_sel[1],"$"), colid ) #exact match
    colid_end <- grep(paste0("^",col_sel[2],"$"), colid)
    stopifnot(colid_start < colid_end)
    col_slice <- colid_start:colid_end
    
    return(list(matSlice = mat[row_slice,col_slice], 
                row_slice=rowid[row_slice],
                col_slice = colid[col_slice]
                )
          )
}
# module1
module1 <- sliceModule(regMat, row_sel = c('Nrip2','C230057M02Rik'), col_sel = c('Gli1','Tcf4'))
module_use = module1
# 
net.dat <- figR.d %>% filter(Score >= 2)
net.dat <- net.dat %>% filter(Motif %in% module_use$col_slice)
net.dat$Motif <- paste0(net.dat$Motif)
net.dat$DORC <- paste0(net.dat$DORC)
tflist.tab <- split(x = net.dat, f = net.dat$Motif)
tflist <- lapply(tflist.tab, FUN = function(x){ x$DORC  } )
# GRN plot
dorcs <- data.frame(name = unique(net.dat$DORC), group = "DORC", size = 3)
tfs <- data.frame(name = unique(net.dat$Motif), group = "TF", size = 8)
nodes <- rbind(dorcs, tfs)
nodes$isinHeatmap <- ifelse(nodes$name %in% module_use$row_slice,'yes','no')
nodes$name <- paste0(nodes$name)
idx <- which(duplicated(nodes$name))
nodes$name[idx] <- paste0( nodes$name[idx],'_',nodes$group[idx]  )
edges <- as.data.frame(net.dat)
cut <- dorcs[grep("Rik$", dorcs$name),]$name
nodes = nodes[!(nodes$name %in% cut),]
edges = edges[!(edges$DORC %in% cut),]
links <- data.frame(
        source = unlist(lapply(edges$Motif, function(x) { which(nodes$name == x) })),
        target = unlist(lapply(edges$DORC, function(x) { which(nodes$name == x) })),
        corr = edges$Corr, 
        enrichment = edges$Enrichment.P
    )
links$Value <- scales::rescale(edges$Score) * 20
links$source_name <- as.character(nodes$name[links$source])
links$target_name <- as.character(nodes$name[links$target])
links <- links[order(links$source_name, links$target_name), ]

g_links <- links[,c('source_name', 'target_name', 'corr')] # TF \ DORCs \ corr
g <- graph_from_data_frame(g_links)
# gene used
vtx <- igraph::V(g)
gene.use <- unique(names(vtx))
n_pathvtx <- length(unique(links$source_name))
# gene expression and DORC score added
dorcmat <- as.data.frame(dorcMat)[gene.use,]
DORCexpre <- data.frame(DORC = log2(rowMeans(dorcmat)+1))
RNAmat = as.data.frame(scale(skin.fib@assays$RNA@data[gene.use,]))
RNAmat[is.na(RNAmat)] <- 0 
RNAexpre <- data.frame(RNA = rowMeans(RNAmat))
expre <- cbind(DORCexpre,RNAexpre)
expre[,"TF"] <- 16
expre[1:n_pathvtx,"TF"] <- 18
# 
color1 = colorRampPalette(colors = c("#5567B1","#B3D5F1","#F8FFEC"))(20)
color2 = colorRampPalette(colors = c("#F8FFEC","#F0686C"))(20) 
x1 <- as.numeric(cut(expre$RNA, breaks = seq(-0.5, 0, length.out = 20)))
x2 <- as.numeric(cut(expre$RNA, breaks = seq(0, 2, length.out = 20)))
palette1 <- color1[x1]
palette2 <- color2[x2]
palette <- c()
for (i in seq_along(palette1)) {
  if (is.na(palette1[i])) {
    palette[i] <- palette2[i]
  } else {
    palette[i] <- palette1[i]
  }
}
V(g)$color <- palette
range(expre$DORC)
map <- function(data, MIN, MAX) {
  d_min <- max(data)
  d_max <- min(data)
  return(MIN + (MAX - MIN) / (d_max - d_min) * (data - d_min))
}
DORCsize <- map(expre$DORC,1,10) * 1.1
V(g)$size <- DORCsize
set.seed(101)
p1= g %>% as_tbl_graph() %>% activate(nodes) %>% mutate(degree  = centrality_degree()) %>% 
      ggraph(layout='nicely') + geom_edge_link(aes(edge_width=links$corr), alpha=.8, colour='grey') + scale_edge_width_continuous(range = c(0.4,1.5)) +
      geom_node_point(aes(size = size), shape=expre$TF, color = cols_f(vcount(g)),show.legend = F) + scale_size(range=c(1, 11) * 1.1) + scale_color_manual(values = cols_f(ecount(g)))+ 
      geom_node_text(aes(label= name), repel=TRUE, size = 3, bg.color = "white") + 
      ggtitle(paste('TF-network ',paste(names(tflist),collapse = '_'  )) ) + theme_void()+ theme(legend.position = 'none') + theme_void()+theme(legend.position = 'none') + theme(legend.position = 'none')+
      coord_equal() 
      #geom_edge_arc(aes(edge_width=links$corr), alpha=.8, colour='grey',strength = 0.05) 
p1