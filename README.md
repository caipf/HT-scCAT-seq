## high throughput single-cell chromatin accessibility and transcriptome sequencing (HT-scCAT-seq)
---
HT-scCAT-seq is a method that simultaneously profiles transcriptome and chromatin accessibility from the same cell.
![[./scCAT-seq.png]]


## Reproducing the results from the paper

```
00_Preprocessing
01_Cellline
02_Brain
03_Skin
04_Skin_Fibroblast
```

## Softwares/Packages needed
**Python packages**
- numpy (v1.19.5)
- scipy (v1.5.3)
- pandas (v1.1.5)
- seaborn (v0.11.0.rc0)
- scanpy (v1.9.0)
- scikit-learn (v0.22.1)

**R packages**
- Matrix (v1.3-3)
- Signac (v1.6.0)
- Seurat (v4.1.0)
- monocle (v2.20.0)
- GenomeInfoDb (v1.30.1)
- EnsDb.Mmusculus.v79 (v2.99.0)
- EnsDb.Hsapiens.v86 (v2.99.0)
- ggplot2 (v3.35)
- patchwork (v1.1.1)
- dplyr (v1.0.8)
- stringr (v1.4.0.9000)
- stringi (v1.7.6)
- FigR (v0.3.5)
- scCustomize (v1.1.3)
- JASPAR2020 (v0.99.10)
- TFBSTools (v1.40.0)
- FNN (v1.1.4.1)
- cisTopic (v2.1.0)
- GenomicRanges (v1.54.1)
- SummarizedExperiment (v1.32.0)
- tidyverse (v2.0.0)
- igraph (2.0.3)
- ComplexHeatmap (v2.18.0)
- monocle3 (v1.3.1)
- cicero (v1.26.0)
- motifmatchr (v1.24.0)
- chromVAR (v1.24.0)

**other Bioinformatics packages**
- chromap (0.2.1-r369)
- MACS2 (v2.2.7.1)
- STARsolo (v2.7.9a)
- Cellranger-ATAC (v2.0.0)
- samtools (v1.9)
- bedtools (v2.30.0)
- seqtk (v1.3-r106)
- sratoolkit (v3.0.0)
- faSize, bedClip, calc, addCols from [UCSC utilities](http://hgdownload.soe.ucsc.edu/admin/exe/)
## Citation
Qiuting D. et al. Deciphering gene regulatory programs in mouse embryonic skin through single-cell multiomics analysis. 
=======
# HT-scCAT-seq
---
