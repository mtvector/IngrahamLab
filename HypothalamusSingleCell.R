## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,cache=TRUE,cache.lazy = FALSE,dev = c('png','pdf'))

## ----CampbelTSNE---------------------------------------------------------
campbell <- FindClusters(object = campbell, reduction.type = "pca", dims.use = 1:20,resolution = .6)
campbell <- StashIdent(object = campbell, save.name = "ClusterNames_0.6")
PCAPlot(campbell,dim.1 = 1, dim.2 = 2,group.by = "ClusterNames_0.6")
#PrintFindClustersParams(object = campbell)
PCAPlot(campbell,dim.1 = 1, dim.2 = 2,group.by="X7.clust_all")
TSNEPlot(campbell,group.by="X7.clust_all")
TSNEPlot(campbell,group.by="ClusterNames_0.6")
FeaturePlot(campbell,features.plot = gnz)
FeaturePlot(campbell,features.plot =gnz)

markers <- FindAllMarkers(
    object = campbell, 
    only.pos = TRUE, 
    min.pct = 0.25, 
    thresh.use = 0.25
)
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

DoHeatmap(
    object = campbell, 
    genes.use = top10$gene, 
    slim.col.label = TRUE, 
    remove.key = TRUE,cex.row=0.1
)

