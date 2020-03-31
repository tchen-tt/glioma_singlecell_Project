# Author chentao
# The function of this script is to draw a diagram
#----------------import pacakges---------------------
library(Seurat)
library(tidyverse)
library(export)
library(RColorBrewer)
library(cowplot)
library(org.Hs.eg.db)
library(clusterProfiler)

#----------------load scRNA-seq data ----------------
load("../ensembl_tumor_health.RData")
load("./cnv.filter.cnvs.RData")

#----------------colors in all pcitures--------------
clust_color <- c("#f8766d", "#de8c00", "#b79f00", "#7cae00",
                 "#00ba38", "#02c086", "#01bfc4", "#00b4f0",
                 "#619cff", "#c77cff", "#f564e3", "#ff64b0")
heatmap_color <- brewer.pal(n = 11, name = "RdBu")
heatmap_color.median <- heatmap_color
heatmap_color.median[6] <- 'lightgray'

#-----------------cell type names ---------------------
Cell_Type <- c("Microglial/Macrophage", "Microglial/Macrophage", 
               "Malignant", "Microglial/Macrophage", "OPCs", "Malignant",
               "Malignant", "Malignant", "Oligodendrocyte", "Astocyte",
               "vascular", "Neuron")
#-----------------define default themes wiht ggplot2---
with.tick.theme <- theme(
  axis.ticks = element_line(color = "black"),
  axis.text = element_text(color = "black", size = 12)
)

#-----------------ploting tSNE-plot--------------------
default.cluster.name <- 0:11  # TSNEplot default cluster name
result.cluster.name <- c(1, 0, 3, 
                         2, 5, 6, 7,
                         4, 8, 9, 10, 11)
#match(default.cluster.name, result.cluster.name)

P.TSNE <- TSNEPlot(ensembl_tumor_health$data) +
  scale_colour_manual(values = clust_color[match(default.cluster.name, result.cluster.name)])
P.TSNE <- P.TSNE + with.tick.theme +
  theme(panel.border = element_rect(color = "black"))

#-----------------ploting marker heatmap---------------
#extract first10 marker gene in every cluster
marker.10.gene <- ensembl_tumor_health$all.marker %>%
  group_by(cluster) %>%
  top_n(., n = 10, wt = avg_logFC)

sample.sort.cluster <- vector() # sample with all cluster
for(i in result.cluster.name) {
  sample.sort.cluster <- c(sample.sort.cluster, WhichCells(ensembl_tumor_health$data, idents = i))
}

genes.sort.cluster <- vector() # top10 marker genes with all cluster
for(i in result.cluster.name) {
  genes.sort.cluster <- c(genes.sort.cluster, 
                          marker.10.gene[marker.10.gene$cluster == i, "gene", drop = TRUE])
}
genes.sort.cluster <- unique(genes.sort.cluster) # make suer genes is unique

# reshape data
# gene cell Expression
P.Heatmap.marker <- ensembl_tumor_health$data[['RNA']]@scale.data[genes.sort.cluster, sample.sort.cluster] %>%
  MinMax(., min = -2.5, max = 2.5) %>%
  as.data.frame %>%
  rownames_to_column('gene') %>%
  gather(., -gene, key = 'cell', value = 'Expression') %>%
  mutate(gene = factor(gene, levels = rev(genes.sort.cluster)),
         cell = factor(cell, levels = sample.sort.cluster)) %>% 
  ggplot(., mapping = aes(x = cell, y = gene, fill = Expression)) +
  geom_raster() + theme(axis.title = element_blank(),
                        axis.text.x = element_blank(),
                        axis.ticks = element_blank()) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradientn(colors = rev(heatmap_color.median))
Heatmap.mark.build <- ggplot_build(P.Heatmap.marker)
Heatmap.mark.y.range <- diff(Heatmap.mark.build$layout$panel_params[[1]]$y.range)
Heatmap.mark.y.pos <- max(Heatmap.mark.build$layout$panel_params[[1]]$y.range) + Heatmap.mark.y.range * 0.015
Heatmap.mark.y.max <- Heatmap.mark.y.pos + Heatmap.mark.y.range * 0.025
cluster.cell.number <- table(ensembl_tumor_health$annotations$tsneClass)
names(cluster.cell.number) <- as.character(0:11)
rester.col <- rep(x = clust_color, times = cluster.cell.number[as.character(result.cluster.name)])

P.Heatmap.marker <- P.Heatmap.marker + annotation_raster(
  raster = t(rester.col), xmin = -Inf, xmax = Inf, ymin = Heatmap.mark.y.pos, ymax = Heatmap.mark.y.max
) + coord_cartesian(ylim = c(0, Heatmap.mark.y.max), clip = "off")

#-----------------------------Scatter plot ---------------------------------------
#marker_gene <- c("ETNPPL", "CALY","SOX10", "CSPG4", "MOG", "IFITM1", "TMEM119",
#                 "AIF1", "EGFR", "TNC")
marker.genes <- c("TMEM119", "AIF1", "EGFR", "TNC", "SOX10", "CSPG4",
                  "MOG", "ETNPPL", "IFITM1", "CALY")
ensembl_tumor_health$annotations <- cbind.data.frame(ensembl_tumor_health$annotations,
                                                     as.data.frame(ensembl_tumor_health$data[['tsne']]@cell.embeddings[rownames(ensembl_tumor_health$annotations),]))
#ensembl_tumor_health$data[['tsne']]@cell.embeddings[rownames(ensembl_tumor_health$annotations),]
P.Scatter <- lapply(marker.genes, FUN = function(x) {
  data <- ensembl_tumor_health$annotations[, c('tSNE_1', 'tSNE_2')]
  data$Expression <- ensembl_tumor_health$data[['RNA']]@data[x, rownames(data), drop = TRUE]
  P.scatter.plot <- data %>% 
    arrange(Expression) %>%
    ggplot(., aes(x = tSNE_1, y = tSNE_2)) +
    geom_point(mapping = aes(color = Expression)) +
    scale_color_gradientn(colors = rev(heatmap_color.median[1:6])) +
    theme_minimal() +
    labs(x = NULL, y = NULL, title = x) +
    theme(plot.title = element_text(hjust = 0.5, face = "italic"),
          axis.ticks = element_blank(),
          axis.text = element_blank())
})
P.Scatter.legend <- suppressMessages(P.Scatter[[4]] + scale_color_gradientn(breaks = c(0, 1, 4),
                                       labels = c('unexpr', 'low', 'high'),
                                       colors = rev(heatmap_color.median[1:6])))
#P.Scatter.legend <- get_legend(P.Scatter.legend)
P.Scatter <- lapply(P.Scatter, FUN = function(x) x + theme(legend.position = "none")) # remove legend
P.Scatter[[length(marker.genes)+1]] <- get_legend(P.Scatter.legend)
plot_grid(plotlist = P.Scatter, nrow = 3, ncol = 4) # plot

#------------------plot immunity variance and tumor cells-----------------------
P.TSNE2 <- ensembl_tumor_health$annotations %>% 
  mutate(., tsneClass = factor(tsneClass, levels = as.character(result.cluster.name))) %>%
  ggplot(., aes(x = tSNE_1, y = tSNE_2)) +
  geom_point(mapping = aes(colour = tsneClass), size = 0.7, shape = 16) + 
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 15),
        panel.background = element_blank(),
        axis.ticks = element_line(color = "black"),
        panel.border = element_rect(fill = NA),
        legend.key = element_blank()) + guides(color = guide_legend(override.aes = list(size = 3))) +
  scale_color_discrete("")
# function get colors
get.color <- function(x = vector()) {
  if(is.null(x)) stop("x is null")
  index <- 1:12
  colors <- clust_color
  colors <- ifelse(test = index %in% x, yes = colors, no = "lightgray")
  return(colors)
}

# plot immunity cluster, tumors culster and vascular cluster
P.TSNE2 + scale_colour_manual(values = get.color(1:3)) + labs(colour = "Class")
P.TSNE2 + scale_colour_manual(values = get.color(4:7)) + labs(colour = "Class")
P.TSNE2 + scale_colour_manual(values = get.color(11))

#---------------------------immunity cluster cells----------------------------
#testing <- FindMarkers(ensembl_tumor_health$data, ident.1 = 1, ident.2 = c(0, 3), min.pct = 0.25, only.pos = TRUE)

# diff sequencet 1, 0, 3 cluster
immu.cel.diff <- lapply(c(1, 0, 3), function(cluster = x, all.cluster = c(1, 0, 3)) {
  FindMarkers(object = ensembl_tumor_health$data, ident.1 = cluster, ident.2 = setdiff(all.cluster, cluster), only.pos = TRUE, min.pct = .25)
})
names(immu.cel.diff) <- as.character(c(1, 0, 3))

#testing.1 <- immu.cel.diff[['3']] %>% 
#  rownames_to_column("gene") %>%
#  arrange(desc(avg_logFC))

#testing <- enrichGO(testing.1$gene[1:50], keyType = "SYMBOL", ont = 'BP', OrgDb = org.Hs.eg.db)
immu.cell.50.gene <- lapply(immu.cel.diff, FUN = function(x) {
  x.sort <- x %>% rownames_to_column('gene') %>% arrange(desc(avg_logFC)) %>%`[[`('gene') %>% 
    `[`(1:50)
#  x.sort <- x.sort[, 'gene', drop = TRUE]
#  x.sort[1:50]
})
#table(testing[[2]] %in% testing[[3]])
#teseting <- ensembl_tumor_health$data[['RNA']]@scale.data[unlist(immu.cell.50.gene),c(WhichCells(ensembl_tumor_health$data, idents = 1),
#                                                                          WhichCells(ensembl_tumor_health$data, idents = 0),
#                                                                          WhichCells(ensembl_tumor_health$data, idents = 3))] %>% 
#  MinMax(., min = -2.5, max = 2.5)
#  heatmap(., Rowv = NA, Colv = NA, scale = "none", col = rev(heatmap_color))
#testing <- ensembl_tumor_health$data[['RNA']]@data[,c(WhichCells(ensembl_tumor_health$data, idents = 1),
#                                                                                     WhichCells(ensembl_tumor_health$data, idents = 0),
#                                                                                     WhichCells(ensembl_tumor_health$data, idents = 3))]

make.diff.genes <- function(cluster = vector(),
                            min.pct = 0.1,
                            only.pos = TRUE,
                            avg_logFc = 0.25) {
  if(is.null(cluster)) stop("input cluster is null")
  if(length(cluster) == 1) top("input more cluster")
  diff.results <- vector(mode = 'list')
  for(i in 1:(length(cluster)-1)) {
    for(j in (i+1):length(cluster)){
      message('Start cluster: ', cluster[i] , ' vs ', cluster[j])
      name <- paste0(cluster[i], '_vs_', cluster[j])
      diff.results[[name]] <- FindMarkers(ensembl_tumor_health$data, ident.1 = cluster[i],
                                          ident.2 = cluster[j],
                                          min.pct = min.pct,
                                          only.pos = only.pos) %>%
        rownames_to_column('gene') %>%
        filter(p_val_adj < 0.05 & avg_logFC >= avg_logFc)
    }
  }
  return(diff.results)
}

make.diff.genes.unique <- function(object = list(),
                                   avg_logFC = 0.25) {
  if(is.null(object)) stop("inputs is unll")
  name <- names(object)
  unique.diff <- vector(mode = 'list', length = length(name))
  for(i in name) {
    other.clus.gene <- vector()
    for(j in setdiff(name, i)) {
      other.clus.gene <- c(other.clus.gene, object[[j]][object[[j]]$avg_logFC >= avg_logFC, 'gene', drop = TRUE])
      other.clus.gene <- unique(other.clus.gene)
    }
    unique.diff[[i]] <- setdiff(object[[i]][object[[i]]$avg_logFC >= avg_logFC, 'gene', drop = TRUE], other.clus.gene)
  }
  return(unique.diff)
}

diff.gene.all <- function(cluster = vector(),
                          logfc.thresholds = 0.25,
                          min.pcts = 0.1,
                          only.pos = TRUE,
                          cut.off.pvale = 0.05) {
  if(is.null(cluster)) top('clusters is null')
  if(length(cluster) == 1) top('input more clusters')
  diff.results <- vector(mode = 'list')
  for(i in cluster) {
    others <- vector(mode = 'list', length = length(cluster) - 1)
    for(j in setdiff(cluster, i)) {
      name <- paste0(i, '_vs_', j)
      message('Start cluster ', name)
      diff <- FindMarkers(ensembl_tumor_health$data, ident.1 = i, ident.2 = j, logfc.threshold = logfc.thresholds,
                  min.pct = min.pcts, only.pos = only.pos)
      diff$gene <- rownames(diff)
      others[[name]] <- diff
    }
    diff.gene <- lapply(others, function(x) {
      x[x$p_val_adj < cut.off.pvale, 'gene', drop = TRUE]
    })
    every.gene <- vector()
    for(k in 1:length(diff.gene)) {
      if(length(every.gene) == 0) every.gene <- diff.gene[[k]]
      every.gene <- intersect(every.gene, diff.gene[[k]])
    }
    diff.results[[as.character(i)]] <- every.gene
  }
  return(diff.results)
}

#----------------deal with immunity cells-------------------------------------
immu.cel.diff <- diff.gene.all(cluster = c(1, 0, 3), min.pcts = 0.6, cut.off.pvale = 0.05, logfc.thresholds = 0.4)
immu.cel.diff.GO <- lapply(immu.cel.diff, FUN = function(x) {
  enrichGO(gene = x, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = 'BP')
})
names(immu.cel.diff.GO) <- paste0('clus', c(1, 0, 3))


immu.cell <- c(WhichCells(ensembl_tumor_health$data, idents = 1),
            WhichCells(ensembl_tumor_health$data, idents = 0),
            WhichCells(ensembl_tumor_health$data, idents = 3))

P.Heatmap.immu <- ensembl_tumor_health$data[['RNA']]@scale.data[unlist(immu.cel.diff), immu.cell] %>% 
  MinMax(., min = -2.5, max = 2.5) %>%
  as.data.frame %>%
  rownames_to_column('gene') %>%
  gather(., -gene, key = 'cell', value = 'Expression') %>%
  mutate(., gene = factor(gene, levels = rev(unlist(immu.cel.diff))),
         cell = factor(cell, levels = immu.cell)) %>%
  ggplot(., aes(x = cell, y = gene)) +
  geom_raster(mapping = aes(fill = Expression)) +
  scale_fill_gradientn(colors = rev(heatmap_color)) +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 4))
heatmap.immu.build <- ggplot_build(P.Heatmap.immu)
heatmap.immu.build.yrange <- diff(heatmap.immu.build$layout$panel_params[[1]]$y.range)
heatmap.immu.y.pos <- max(heatmap.immu.build$layout$panel_params[[1]]$y.range) + heatmap.immu.build.yrange * 0.015
heatmap.immu.y.max <- heatmap.immu.y.pos + heatmap.immu.build.yrange * 0.025
P.Heatmap.immu <- P.Heatmap.immu + annotation_raster(raster = t(rep(clust_color[1:3], times = table(ensembl_tumor_health$annotations$tsneClass)[1:3])),
                                   xmin = -Inf, xmax = Inf,
                                   ymin = heatmap.immu.y.pos, ymax = heatmap.immu.y.max) +
  coord_cartesian(ylim = c(0, heatmap.immu.y.max), clip = "off")

#--------------go terms --------------------------------
GOID <- list(clus1 = c("GO:0050867", "GO:0045785", "GO:0042063",
                       "GO:0014015", "GO:1905517", "GO:0031663",
                       "GO:0060326", "GO:1903975", "GO:0042089",
                       "GO:1901214"),
             clus0 = c("GO:0019882", "GO:0060333", "GO:0034341", 
                       "GO:0050851", "GO:0002694", "GO:0050867", 
                       "GO:0060192"),
             clus3 = c("GO:0045047", "GO:0061621", "GO:0006007", 
                       "GO:0006734", "GO:0042119", "GO:0009435", 
                       "GO:0019320", "GO:0046365", "GO:0006096",
                       "GO:0046031", "GO:0071456", "GO:0001666"))
# modify clus1, clus0, clus3
P.immue.GO <- lapply(X = c(1, 0, 3), function(x) {
  index <- x
  index <- ifelse(x == 0, yes = 2, no = x)
  p.immue.GO <- immu.cel.diff.GO[[paste0('clus', x)]]@result %>%
    filter(ID %in% GOID[[paste0('clus', x)]]) %>% 
    ggplot(., aes(x = reorder(Description, -log10(p.adjust)), y = -log10(p.adjust))) +
    geom_bar(stat = 'identity', fill = clust_color[index], alpha = 1) + coord_flip() +
    theme(panel.background = element_blank(),
          axis.line.x = element_line(),
          axis.ticks.y = element_blank(),
          axis.ticks = element_line(color = "black"),
          axis.text.x = element_text(size = 12, color = 'black'),
          axis.text.y = element_text(size = 12, color = 'black'),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 15, color = 'black')) +
    scale_y_continuous(expand = c(0, 0))
  return(p.immue.GO)
})
P.IMMUE.GO.COB <- plot_grid(plotlist = P.immue.GO, nrow = 3, rel_heights =c(10, 7, 12), align = "v")


#Normal.cell <- diff.gene.all(cluster = c(0, 4, 8, 9, 11), min.pcts = 0.4, only.pos = TRUE, logfc.thresholds = 0.3)
#Cancer.cell <- diff.gene.all(cluster = c(2, 5, 6))



#--------------------------------Let's go Cancer cells----------------------------------------
Cancer.cell <- WhichCells(ensembl_tumor_health$data, idents = c(2, 5, 6, 7))
Cancer.cell.data <- ensembl_tumor_health$data[['RNA']]@data[, Cancer.cell]
Cancer.cell.diff <- diff.gene.all(cluster = c(2, 5, 6, 7), logfc.thresholds = 0.3, min.pcts = 0.5)

MES.like.gene <- c("DDIT3", "ENO2", "VIM", "ADM", "LDHA", "HILPDA", "ANXA1", "ANXA2",
                   "CHI3L1", "CD44", "RELB", "TRADD", "PDPN")
NPC.like.gene <- c("SOX4", "DCX", "DLL3", "SOX11", "RND3", "STMN4", "STMN2",
                   "DLX5", "DLX6", "NCAM1", "NKX2-5", "ASCL1")
ACS.like.gene <- c("CST3", "GFAP", "S100B", "HOPX", "SLC1A3", "MLC1", "AQP4",
                   "SPARCL1")
OPC.like.gene <- c("PLP1", "ALCAM", "OLIG1", "OMG", "PLLP",
                   "OLIG2", "TNR", "SMOC1", "PDGFRA", "BCAN", "CCND1",
                   "CSPG4")


Four.like.gene <- list(MES = MES.like.gene,
                       NPC = NPC.like.gene,
                       ACS = ACS.like.gene,
                       OPC = OPC.like.gene)

lapply(Cancer.cell.diff, function(x) table(x %in% NPC.like.gene))

Cancer.testing <- FindMarkers(ensembl_tumor_health$data,
                              ident.1 = 7, ident.2 = c(2, 5, 6),
                              min.pct = 0.5, logfc.threshold = 0.3) %>%
  rownames_to_column('gene')
#Cancer.testing.1 <-  Cancer.testing %>% top_n(n = 50, wt = avg_logFC)
#ts.go <- enrichGO(gene = Cancer.cell.diff[['2']], OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP")

cencer.diff <- diff.gene.all(cluster = c(2, 5, 6),
                             min.pcts = 0.4,
                             logfc.thresholds = 0.3)

names(cencer.diff)
#ts.go <- enrichGO(gene = cencer.diff[['6']], OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP")

sapply(X = Four.like.gene, FUN = function(clus) {
  sapply(cencer.diff, FUN = function(x) length(intersect(x, clus)))
})

make.diff.cluster <- function(object, 
                              cluster = vector(),
                              min.pct = 0.4,
                              logfc.thred = 0.3,
                              p.adj = 0.05) {
  diff.result <- vector(mode = "list")
  for(i in cluster) {
    message("Star analysis cluster ", i)
    cluster.other <- setdiff(cluster, i)
    diff.result[[as.character(i)]] <- FindMarkers(object = object, 
                                                  ident.1 = i, ident.2 = cluster.other,
                                                  min.pct = min.pct, logfc.threshold = logfc.thred,
                                                  only.pos = TRUE) %>%
      rownames_to_column('gene') %>%
      filter(p_val_adj < p.adj) %>%
      arrange(desc(avg_logFC))
  }
  return(diff.result)
}

testing <- make.diff.cluster(object = ensembl_tumor_health$data, 
                             cluster = c(2, 5, 6, 7),
                             min.pct = 0.4,
                             logfc.thred = 0.3)
sapply(X = Four.like.gene, FUN = function(clus) {
  sapply(testing, FUN = function(x) length(intersect(x$gene, clus)))
})

ts.go <- enrichGO(gene = testing[["6"]]$gene, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP")
kegg.gene <- mapIds(x = org.Hs.eg.db, keys = testing[["6"]]$gene, column = "UNIPROT", keytype = "SYMBOL", multiVals = "first")
ts.kegg <- enrichMKEGG(gene = kegg.gene, organism = "hsa", keyType = "uniprot")




t <- sapply(Four.like.gene, FUN = function(x) {
  data <- ensembl_tumor_health$data[['RNA']]@data[x, WhichCells(ensembl_tumor_health$data, idents = c(2, 5, 6, 7))]
  data <- t(scale(t(data)))
  resunts <- apply(data, 2, FUN = function(x) {
    mean(sort(x, decreasing = TRUE)[1:5])
  })
})

tumor.sample.tsne.class <- ensembl_tumor_health$annotations %>%
  rownames_to_column('sample') %>%
  dplyr::filter(tsneClass %in% c(2, 5, 6, 7))

t %>% as.data.frame %>%
  rownames_to_column('sample') %>%
  gather(-sample, key = 'type', value = 'expression') %>%
  left_join(., tumor.sample.tsne.class, by = 'sample') %>%
  group_by(., type, tsneClass) %>%
  summarise(type.mean = mean(expression),
            type.sem = sd(expression)/sqrt(n()),
            n = n()) %>%
  ggplot(., aes(x = tsneClass, y = type.mean, fill = type)) +
  geom_bar(stat = 'identity', position = position_dodge(width=0.9)) + theme_classic() +
  geom_errorbar(aes(ymin = type.mean - type.sem, ymax = type.mean + type.sem), position = position_dodge(width=0.9), width = 0.2) +
  labs(y = "Relative Expression") + scale_y_continuous(expand = c(0, 0.05)) +
  scale_fill_manual(values = unique(GENE$gene_color))
ggsave("single.cell.fourtype.classificy.pdf", width = 6, height = 4)


annotations <- ensembl_tumor_health$annotations %>%
  rownames_to_column('sample') %>%
  filter(sample %in% filter.cnvs$smaple)

ensembl_tumor_health$data[['RNA']]@data[c("HIF1A", 'EPAS1'), annotations$sample] %>% 
  as.data.frame %>% rownames_to_column('gene') %>%
  gather(-gene, key = 'sample', value = 'expression') %>%
  left_join(., annotations) %>% 
  mutate(tsneClass = factor(tsneClass, levels = c(2, 7, 5, 6))) %>%
  filter(gene != 'HIF1A' & expression > 0) %>%  #EPAS1
  ggplot(., aes(x = tsneClass, y = expression)) + geom_boxplot(aes(fill = tsneClass)) + theme_classic() +
  geom_jitter() + theme(legend.position = 'none') + labs(title = "HIF2A") +
  theme(plot.title = element_text(hjust = 0.5, face = "italic"))
ggsave(filename = "HIF2A.plot.pdf", width = 3.25, height = 2.7)









as.1 <- ensembl_tumor_health$all.marker %>% group_by(cluster) %>% filter(cluster == 9) %>% top_n(n = 100, wt = avg_logFC) %>% .[, "gene", drop = TRUE]
as <- FindMarkers(ensembl_tumor_health$data, ident.1 = 11, ident.2 = c(9, 10, 4),
                  min.pct = 0.5, logfc.threshold = 0.3)
as.2 <- as %>% rownames_to_column(., 'gene') %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_logFC)) %>%
  .[1:100, 'gene']
tts <- GENE %>%
  filter(Subtype == "NEURALSUBTYPE") %>% .[,'SYMBOL']
DoHeatmap(object = ensembl_tumor_health$data, features = as.2)
#
date()




testing <- FindMarkers(ensembl_tumor_health$data, ident.1 = 6, ident.2 = c(2, 5, 7),
                       min.pct = 0.4, logfc.threshold = 0.3)
t.6 <- testing %>%
  rownames_to_column('gene') %>% 
  arrange(desc(avg_logFC)) %>%
  top_n(n = 150, wt = -avg_logFC)
testing.5 <- FindMarkers(ensembl_tumor_health$data, ident.1 = 5, ident.2 = 6,
                         min.pct = 0.4, logfc.threshold = 0.3)

t.5 <- testing.5 %>%
  rownames_to_column('gene') %>% 
  arrange(desc(avg_logFC)) %>%
  top_n(n = 150, wt = avg_logFC)




T.testing <- FindMarkers(ensembl_tumor_health$data,
                         ident.1 = 2, ident.2 = 7, min.pct = 0.4, logfc.threshold = 0.3)

T.testing %<>% rownames_to_column('gene') %>%
  arrange(desc(avg_logFC)) %>% filter(p_val_adj < 0.05)
TT <- T.testing[T.testing$avg_logFC>0,]$gene
m <- enrichGO(gene = TT, keyType = "SYMBOL", ont = "BP", OrgDb = org.Hs.eg.db)


TT.L <- T.testing[T.testing$avg_logFC<0,]$gene
m.L <- enrichGO(gene = TT.L, keyType = "SYMBOL", ont = "BP", OrgDb = org.Hs.eg.db)


#
get.colrs.alpha <- function(colors, alpha = 0.7) {
  colors <- t(col2rgb(colors) * alpha / 255)
  color <- rgb(red = colors[, 1], green = colors[, 2], blue = colors[, 3], maxColorValue = 1)
  return(color)
}


colors <- function (colour, alpha = NA) 
{
  col <- grDevices::col2rgb(colour, TRUE)/255
  if (length(colour) != length(alpha)) {
    if (length(colour) > 1 && length(alpha) > 1) {
      stop("Only one of colour and alpha can be vectorised")
    }
    if (length(colour) > 1) {
      alpha <- rep(alpha, length.out = length(colour))
    }
    else if (length(alpha) > 1) {
      col <- col[, rep(1, length(alpha)), drop = FALSE]
    }
  }
  alpha[is.na(alpha)] <- col[4, ][is.na(alpha)]
  new_col <- grDevices::rgb(col[1, ], col[2, ], col[3, ], 
                            alpha)
  new_col[is.na(colour)] <- NA
  new_col
}






#
P <-TSNEPlot(ensembl_tumor_health$data) + 
  theme(panel.border = element_rect(color = "black", size = 1),
        axis.text = element_text(size = 12, color = "black"),
        axis.line = element_blank())
graph2ppt(file = "Figure1.A.pptx")

###########annotaitons##################################
#
#cluster 6 PN, Cluster 2 MES, Cluster 5 Astocyte
#cluster 7 OPCS like
#########################################################
Cell_Type <- c("Microglial/Macrophage", "Microglial/Macrophage", 
               "Malignant", "Microglial/Macrophage", "OPCs", "Malignant",
               "Malignant", "Malignant", "Oligodendrocyte", "Astocyte",
               "vascular", "Neuron")
names(Cell_Type) <- 0:11
ensembl_tumor_health$annotations$tsneClass <- as.numeric(ensembl_tumor_health$annotations$tsneClass)
ensembl_tumor_health$annotations$Cell_Type <- Cell_Type[ensembl_tumor_health$annotations$tsneClass]

cell_number <- table(ensembl_tumor_health$annotations$Cell_Type, ensembl_tumor_health$annotations$tissue) %>%
  as.data.frame %>%
  spread(., key = Var2, value = Freq) %>%
  rename(., `Cell Type` = Var1) %>%
  mutate(`Total of Cells` = Periphery + Tumor) %>%
  mutate(Percentage = round(`Total of Cells`/ sum(`Total of Cells`) * 100))
cell_number$Percentage <- paste0(cell_number$Percentage, "%")  
colnames(cell_number) <- c("Cell Type", "Periphery", "Tumor", "Total of Cells", "Percentage")
  
top5_gene <- ensembl_tumor_health$all.marker %>% group_by(cluster) %>%
  top_n(n = 5, wt = avg_logFC)
all.marker <- FindAllMarkers(ensembl_tumor_health$data, min.pct = 0.25)
top20_gene <- ensembl_tumor_health$all.marker %>% group_by(cluster) %>%
  top_n(n = 20, wt = avg_logFC)
top10_gene <- all.marker %>% 
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_logFC)
top50_gene <- ensembl_tumor_health$all.marker %>% group_by(cluster) %>%
  top_n(n = 50, wt = avg_logFC)
t.3 <- FindMarkers(ensembl_tumor_health$data, ident.1 = 3, ident.2 = c(1, 0),
            only.pos = TRUE, min.pct = 0.2) %>%
  rownames_to_column("gene") %>% top_n(n = 50, wt = avg_logFC)
t.1 <- FindMarkers(ensembl_tumor_health$data, ident.1 = 1, ident.2 = c(3, 0))

plots <- lapply(c(0, 1, 3, 2, 5, 6, 7, 4, 8, 9, 10, 11), function(x){
  DoHeatmap(ensembl_tumor_health$data, features = top5_gene$gene,
            cells = WhichCells(ensembl_tumor_health$data, idents = x)) + NoLegend() + theme(axis.text.y = element_blank())
})
P <- DoHeatmap(ensembl_tumor_health$data, features = top5_gene$gene, cells = 1000)
sort_level <- factor(Idents(ensembl_tumor_health$data), levels = c(0, 1, 3, 2, 5, 6, 7, 4, 8, 9, 10, 11)) %>% sort()
P$data$Cell <- factor(P$data$Cell, levels = names(sort_level))



plots <- mtcars %>%
  rownames_to_column("company") %>%
  gather(-company, key = "charactors", value = "level") %>%
  ggplot(., aes(x = company, y = charactors, fill = level)) +
  geom_raster() + theme(axis.text = element_blank(),
                      axis.ticks = element_blank(),
                      axis.title = element_blank())
plots.buld <- ggplot_build(plots)
plots.buld$layout$panel_params[[1]]$y.range
y.range <- diff(x = plots.buld$layout$panel_params[[1]]$y.range)
y.pos <- max(plots.buld$layout$panel_params[[1]]$y.range) + y.range * 0.015
y.max <- y.pos + 0.025 * y.range
P1 <- plots + annotation_raster(raster = "red", xmin = -Inf, xmax = Inf, ymin = y.pos,
                          ymax = y.max) +
  coord_cartesian(ylim = c(0, y.max), clip = "off") + guides(fill = "none") +
  scale_x_discrete(expand = c(0, 0)) + theme(panel.background = element_blank())

white.rect = element_rect(fill = 'red')


P1 + geom_rect(., aes(x = 1, y = 3))




rm(top5_gene)

annotaiton.sample <- ensembl_tumor_health$annotations %>% 
  rownames_to_column("sample") %>%
  mutate(., tsneClass = factor(tsneClass, 
                               levels = c(1, 0, 3, 2, 6, 5, 7, 4, 8, 9, 10, 11))) %>%
  arrange(tsneClass)
sample.genes <- top10_gene %>%
  ungroup(cluster) %>%
  mutate(., cluster = factor(cluster, 
                             levels = c(1, 0, 3, 2, 6, 5, 7, 4, 8, 9, 10, 11))) %>%
  arrange(cluster)
P.top10.gene <- ensembl_tumor_health$data[['RNA']]@scale.data[unique(sample.genes$gene), annotaiton.sample$sample] %>%
  MinMax(., min = -2.5, 2.5) %>%
  as.data.frame %>%
  rownames_to_column("gene") %>%
  gather(-gene, key = 'cell', value = 'Expression') %>%
  mutate(gene = factor(gene, levels = rev(unique(sample.genes$gene))),
         cell = factor(cell, levels = annotaiton.sample$sample)) %>%
  ggplot(., aes(x = cell, y = gene)) +
  geom_raster(mapping = aes(fill = Expression))
P.top10.gene <- P.top10.gene + theme(axis.title = element_blank(),
                     axis.ticks = element_blank(),
                     axis.text.x = element_blank()) +
  scale_fill_gradientn(colors = rev(brewer.pal(11, "RdBu"))) +
  theme(legend.position = "bottom",
        axis.text.y = element_blank())


GBM <- FindMarkers(ensembl_tumor_health$data, ident.1 = c(2, 5, 6, 7),
                   ident.2 = c(0, 1, 3, 4, 8, 9, 10, 11))
GBM %>% 
  rownames_to_column("gene") %>%
  arrange(desc(avg_logFC))


#-------------------------immunity cluster ------------------------------------

immunity_color <- clust_color
immunity_color <- ifelse(test = immunity_color %in% c("#f8766d", "#de8c00", "#7cae00"),
                         yes = immunity_color,
                         no = "gray")
TSNEPlot(ensembl_tumor_health$data) + 
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black"),
        panel.border = element_rect()) + 
  scale_colour_manual(values = immunity_color)

P.immune <- TSNEPlot(ensembl_tumor_health$data) + 
  theme(panel.border = element_rect(color = "black", size = 1),
        axis.text = element_text(size = 12, color = "black"),
        axis.line = element_blank()) + scale_colour_manual(values = immunity_color)
graph2pdf(x = P.immune, file = "Figure2.immunoity.pdf")

heatmap_color <- brewer.pal(11, "RdBu")
heatmap_color[6] <- "gray"
DoHeatmap(ensembl_tumor_health$data, features = top10_gene$gene, draw.lines = FALSE, cells = WhichCells(ensembl_tumor_health$data, idents = 0)[1:240]) +
  scale_fill_gradientn(colours = rev(brewer.pal(11, "RdBu")))

#------------------do diff expression in tumor cells ------------------------------------
do.diff.expression <- function(objects,
                               cell.cluster1 = vector(),
                               cell.cluster2 = vector(),
                               min.pct = 0.4,
                               logfc.threshold = 0.3,
                               only.pos = TRUE,
                               ...) {
  diff.results <- FindMarkers(object = objects, 
                              cells.1 = cell.cluster1, cells.2 = cell.cluster2, 
                              min.pct = min.pct, logfc.threshold = logfc.threshold, only.pos = only.pos)
  diff.results %<>% rownames_to_column(., 'gene')
  return(diff.results)
}

annotations.2_7 <- ensembl_tumor_health$annotations %>%
  rownames_to_column('smaple') %>%
  filter(tsneClass %in% c(2, 5, 6, 7) & smaple %in% filter.cnvs$smaple)

filter.cnvs %<>%
  left_join(., annotations.2_7, by = 'smaple')

class.2 <- do.diff.expression(ensembl_tumor_health$data[['RNA']]@data,
                              cell.cluster1 = filter.cnvs$smaple[filter.cnvs$tsneClass == 2],
                              cell.cluster2 = filter.cnvs$smaple[filter.cnvs$tsneClass %in% c(5, 6)])
class.7 <- do.diff.expression(ensembl_tumor_health$data[['RNA']]@data,
                              cell.cluster1 = filter.cnvs$smaple[filter.cnvs$tsneClass == 7],
                              cell.cluster2 = filter.cnvs$smaple[filter.cnvs$tsneClass %in% c(5, 6)])
class.5 <- do.diff.expression(ensembl_tumor_health$data[['RNA']]@data,
                              cell.cluster1 = filter.cnvs$smaple[filter.cnvs$tsneClass == 5],
                              cell.cluster2 = filter.cnvs$smaple[filter.cnvs$tsneClass %in% c(2, 6, 7)])
class.6 <- do.diff.expression(ensembl_tumor_health$data[['RNA']]@data,
                              cell.cluster1 = filter.cnvs$smaple[filter.cnvs$tsneClass == 6],
                              cell.cluster2 = filter.cnvs$smaple[filter.cnvs$tsneClass %in% c(5, 7, 2)])
class.2.50 <- class.2 %>% filter(p_val_adj < 0.05) %>%
  top_n(n = 50, wt = avg_logFC)
class.7.50 <- class.7 %>% filter(p_val_adj < 0.05) %>%
  top_n(n = 50, wt = avg_logFC)
class.5.50 <- class.5 %>% filter(p_val_adj < 0.05) %>%
  top_n(n = 50, wt = avg_logFC)
class.6.50 <- class.6 %>% filter(p_val_adj < 0.05) %>%
  top_n(n = 50, wt = avg_logFC)

filters.gene <- c()

P.tumors.heatmap <- function(genes = vector(),
                             samples = vector(),
                             ...) {
  ensembl_tumor_health$data[['RNA']]@scale.data[genes, samples] %>%
    MinMax(min = -2.5, max = 2.5) %>%
    as.data.frame %>%
    rownames_to_column('gene') %>%
    gather(-gene, key = "sample", value = "expression") %>%
    mutate(sample = factor(sample, levels = samples),
           gene = factor(gene, levels = rev(genes))) %>%
    ggplot(., aes(x = sample, y = gene)) +
    geom_raster(mapping = aes(fill = expression)) +
    scale_fill_gradientn(colors = rev(heatmap_color)) +
    theme(axis.ticks = element_blank(),
          axis.text = element_blank(), legend.position = "bottom")
    
    
}

tumors.filter.orders <- vector()
for(i in c(2, 7, 5, 6)) {
  tumors.filter.orders <- c(tumors.filter.orders, filter.cnvs$smaple[filter.cnvs$tsneClass == i])
}
P.tumors.heatmap(genes = c(class.7.50$gene,class.5.50$gene, class.6.50$gene), samples = tumors.filter.orders)
ggsave(filename = "class7.6.5.pdf", width = 5.5, height = 9.5)
P.tumors.heatmap(genes = c(class.2.50$gene,class.5.50$gene, class.6.50$gene), samples = tumors.filter.orders)
ggsave(filename = "class2.6.5.pdf", width = 5.5, height = 9.5)



