#paths
args <- commandArgs()
baseName <- args[6]

library(Seurat,lib = "/Users/lihuazhang/Documents/Rlibrary")# path to Seurat v2
library(dplyr)
library(plyr)
library(data.table)
#sessionInfo()

# combine and preprocessing
infile <- paste(baseName, "data_temporal.txt", sep="/")
data <- read.table(file = infile,sep = '\t',row.names=1,header=T)
w10x_new <- CreateSeuratObject(raw.data = data, min.cells = 3, min.genes = 200, project = "MUT")

# Normalizing the data
w10x_new <- NormalizeData(object = w10x_new, normalization.method = "LogNormalize", scale.factor = 10000)

# Detection of variable genes
w10x_new <- FindVariableGenes(object = w10x_new, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.1, x.high.cutoff = 5, y.cutoff = 0.25)
length(w10x_new@var.genes)

w10x_new <- ScaleData(object = w10x_new)

### Perform linear dimensional reduction
w10x_new <- RunPCA(w10x_new, pc.genes = w10x_new@var.genes, pcs.compute = 30, do.print = FALSE)
#PCAPlot(object = w10x_new, dim.1 = 1, dim.2 = 2)

numPC <- 30
res <- c(0.1,0.2,0.3,0.4,0.5)
n <- length(res)
for (i in 1:n)
{
  w10x_new <- FindClusters(object = w10x_new, reduction.type = "pca", dims.use = 1:numPC, resolution = res[[i]],algorithm = 1,save.SNN = TRUE,print.output = 0,force.recalc = T)
  print(length(unique(w10x_new@ident)))
  out <- paste(baseName, paste("identity_clustering_res",res[[i]],".txt",sep = ""), sep="/")
  write.table(w10x_new@ident,file = out,sep = '\t')
}



