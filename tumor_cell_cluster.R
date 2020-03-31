# Created on: 2019-12-14
tumorsample <- WhichCells(gliomass, idents = c(2, 5, 6, 7))
Peripherysample <- annotations %>% rownames_to_column("samples") %>% 
  filter(tissue == "Periphery" & tsneClass %in% c(2, 5, 6, 7))
corsamples <- setdiff(tumorsample, Peripherysample$samples)


CancerCounts <- gliomass[["RNA"]]@counts
CancerCounts <- CreateSeuratObject(counts = CancerCounts, project = "cancer")
CancerCounts <- NormalizeData(CancerCounts, normalization.method = "LogNormalize", scale.factor = 10000)

CancerCounts <- FindVariableFeatures(CancerCounts, selection.method = "vst", nfeatures = 3000)
Cnacertop20 <- head(VariableFeatures(CancerCounts), 20)
plotA <- VariableFeaturePlot(CancerCounts)
plotB <- LabelPoints(plot = plotA, points = Cnacertop20)
CombinePlots(plots = list(plotA, plotB), ncol = 2, legend = "bottom")

Cancerall.genes <- rownames(CancerCounts)
CancerCounts <- ScaleData(CancerCounts, features = Cancerall.genes)

CancerCounts <- RunPCA(CancerCounts, features = VariableFeatures(object = CancerCounts))
DimPlot(CancerCounts, reduction = "pca")
DimHeatmap(CancerCounts, dims = 1, cells = ncol(CancerCounts), balanced = TRUE)

CancerCounts <- JackStraw(CancerCounts, num.replicate = 100)
CancerCounts <- ScoreJackStraw(CancerCounts, dims = 1:20)
JackStrawPlot(CancerCounts, dim = 1:15)
ElbowPlot(CancerCounts, ndims = 100)

CancerCounts <- FindNeighbors(CancerCounts, dims = 1:20)
CancerCounts <- FindClusters(CancerCounts, resolution = 0.08)

CancerAllMarker <- FindAllMarkers(CancerCounts)
CancerAllMarker %<>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
DoHeatmap(CancerCounts, features = CancerAllMarker$gene)
CancerCounts <- RunTSNE(CancerCounts, dims = 1:10)

genes <- mapIds(org.Hs.eg.db, keys = GENE$ENSEMBL, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first") %>% as.data.frame %>% set_colnames("gene") %>% rownames_to_column("ENSEMBL")
genes$gene <- as.character(genes$gene)
DoHeatmap(CancerCounts, features = genes$gene)
DoHeatmap(CancerCounts, features = c("DDIT3", "ENO2", "VIM", "ADM", "LDHA", "HILPDA",
                                     "VIM", "ANXA1", "ANXA2", "CHI3L1", "CD44",
                                     "CST3", "GFAP", "S100B", "HOPX", "SLC1A3", "MLC1",
                                     "PLP1", "ALCAM", "OLIG1", "OMG", "PLLP"))

gliomass_scale <- gliomass[["RNA"]]@scale.data
gliomass_scale[gliomass_scale > 2] <- 2
gliomass_scale[gliomass_scale < -2] <- -2
gliomass_scale[, Peripherysample$samples]
genes$gene
heatmap(gliomass_scale[genes[genes$gene %in% rownames(gliomass_scale),]$gene, Peripherysample$samples],
        col = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
        Rowv = NA, scale = "none", margins = c(10, 4))

pheatmap::pheatmap(gliomass_scale[genes[genes$gene %in% rownames(gliomass_scale),]$gene, Peripherysample$samples],
                   color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
                   cluster_rows = TRUE)


mark.5 <- FindMarkers(CancerCounts, ident.1 = 5, ident.2 = c(1, 4))
mark.5 %<>% rownames_to_column("gene")
t.5 <- mark.5 %>% top_n(n = 150, wt = avg_logFC)

mark.1 <- FindMarkers(CancerCounts, ident.1 = 1, ident.2 = c(4, 5))
mark.1 %<>% rownames_to_column("gene")
t.1 <- mark.1 %<>% top_n(n = 150, wt = avg_logFC)

makr.4 <- FindMarkers(CancerCounts, ident.1 = 4, ident.2 = c(1, 5))
makr.4 %<>% rownames_to_column("gene")
t.4 <- makr.4 %>% top_n(n = 150, wt = avg_logFC)


all.marker <- FindAllMarkers(CancerCounts)
maks <- all.marker %>% group_by(cluster) %>%
  top_n(n = 100, wt = avg_logFC)
table(t.4$gene %in% mt$gene)
