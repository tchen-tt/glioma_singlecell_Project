temporary <- list()
temporary$expr_count <- counts
temporary$infer_gene <- gene_order_file
temporary$infer_gene <- read.table(temporary$infer_gene, sep = "\t", stringsAsFactors = FALSE)
temporary$gtf <- gtf
temporary$gliomas_untnse <- gliomas
temporary$gliomas_tnse <- gliomass
temporary$all_marker <- gliomass.markers
temporary$sample_annotation <- annotations

get.samples <- function(x = vector()) {
  sampless <- vector()
  clusters <- vector()
  index <- 12
  for(i in x) {
    extract.data <- annotations[annotations$tsneClass == i, ]
    dieta <- ifelse(nrow(extract.data) >= 100, yes = 0.1, no = 0.5)
    set.seed(10001)
    sampless <- c(sampless, sample(x = rownames(extract.data), replace = FALSE, size = ceiling(nrow(extract.data) * dieta)))
    clusters <- c(clusters, rep(as.character(index), ceiling(nrow(extract.data) * dieta)))
    index <- index + 1
  }
  names(clusters) <- sampless
  return(clusters)
}
annotations.copy <- annotations
random.sample <- get.samples(x = as.character(x = c(1, 0, 3, 4, 8, 9, 10, 11)))
annotations.copy$tsneClass <- as.character(annotations.copy$tsneClass)
annotations.copy[names(random.sample),]$tsneClass <- random.sample


tt <- GeteCnv(temporary$gliomas_tnse,
              gtfs = gtf$uniquegene[gtf$uniquegene$GENENAMES %in% temporary$infer_gene$V1,],
              annotations = annotations.copy,
              chr = c("Y", "MT"),
              refgroup = paste0('class', 12:19),
              outputdir = "/home/taotao/scRNA/project/shiny/scRNA_gliomas/cnv_gliomass1220.every.normal",
              scale.data = FALSE)


image(log2(ll@expr.data)[, WhichCells(gliomass, idents = 2)], 
      col = rev(brewer.pal(11, "RdBu")), 
      breaks = c(-0.884,-0.78, seq(-0.58,-0.09, length.out = 4), seq(0.09,0.55, length.out = 4),0.7, 0.988))


index <- sample(1:nrow(annotations), size = nrow(annotations), replace = FALSE)
annotations <- annotations[index,]

image(log2(ll@expr.data)[,GetClustSample(annotations, cluster = c(4, 8:11))], 
      col = rev(brewer.pal(11, "RdBu")), 
      breaks = c(-0.884,-0.78, seq(-0.58,-0.09, length.out = 4), seq(0.09,0.55, length.out = 4),0.7, 0.988),
      axes = F)

image(log2(ll@expr.data)[,GetClustSample(annotations, cluster = c(1, 0, 3))], 
      col = rev(brewer.pal(11, "RdBu")), 
      breaks = c(-0.884,-0.78, seq(-0.58,-0.09, length.out = 4), seq(0.09,0.55, length.out = 4),0.7, 0.988),
      axes = F)

image(log2(ll@expr.data)[,GetClustSample(annotations, cluster = c(2, 5:7))], 
      col = rev(brewer.pal(11, "RdBu")), 
      breaks = c(-0.884,-0.78, seq(-0.58,-0.09, length.out = 4), seq(0.09,0.55, length.out = 4),0.7, 0.988),
      axes = F)


vhline <- function( h =  NULL) {
  if(!is.null(h)) {
    print("points lines")
  }
  abline(v = a, lwd = 0.7, col = "gray80")
}





breaks = c(-0.884,-0.7, seq(-0.58,-0.09, length.out = 4), seq(0.09,0.55, length.out = 4),0.6, 0.988)

opar <- par(no.readonly = TRUE)
layout(matrix(1:3, nrow = 3), heights = c(2, 6, 5))
par(mai = c(0.1, 1.5, 0.1, 1.5))
image(1:13618, 1:621,log2(ll@expr.data)[,GetClustSample(annotations, cluster = c(4, 8:11))], 
      col = rev(brewer.pal(11, "RdBu")), 
      breaks = c(-0.884,-0.78, seq(-0.58,-0.09, length.out = 4), seq(0.09,0.55, length.out = 4),0.7, 0.988),
      axes = F, ann = F)

par(mai = c(0.1, 1.5, 0.1, 1.5))
image(1:13618, 1:1842,log2(ll@expr.data)[,GetClustSample(annotations, cluster = c(1, 0, 3))], 
      col = rev(brewer.pal(11, "RdBu")), 
      breaks = c(-0.884,-0.78, seq(-0.58,-0.09, length.out = 4), seq(0.09,0.55, length.out = 4),0.7, 0.988),
      axes = F, ann = F)

par(mai = c(1, 1.5, 0.1, 1.5))
image(1:13618, 1:1103, log2(ll@expr.data)[,GetClustSample(annotations, cluster = c(2, 5:7))], 
      col = rev(brewer.pal(11, "RdBu")), 
      breaks = c(-0.884,-0.78, seq(-0.58,-0.09, length.out = 4), seq(0.09,0.55, length.out = 4),0.7, 0.988),
      axes = F, ann = F)
par(opar)


tiff(filename = "./cnv_heatmap.tif", width = 1111*5, height = 739*5, res = 72*5)
opar <- par(no.readonly = TRUE)
layout(matrix(1:3, nrow = 3), heights = c(2, 6, 5))
par(mai = c(0.1, 2, 0.1, 1.5))
image(1:13618, 1:621,log2(ll@expr.data)[,GetClustSample(annotations, cluster = c(4, 8:11))], 
      col = rev(brewer.pal(11, "RdBu")), 
      breaks = c(-0.884, -0.6, -0.4, -0.25, -0.2, -0.1, 0.1, 0.3, 0.4, 0.45, 0.6, 0.988),
      axes = F, ann = F)
abline(v = cumsum(table(ll@gene_order$chr))[-23], col = "gray60", lwd = 0.8)
abline(h = cumsum(table(annotations$tsneClass)[c(4, 8:11)+1])[-5], col = "gray60", lwd = 0.8)
axis(1, col = "white", at = cumsum(table(ll@gene_order$chr))[-23], labels = NA, col.ticks = "black", tck = 0.04)
axis(3, col = "white", at = cumsum(table(ll@gene_order$chr))[-23], labels = NA, col.ticks = "black", tck = 0.04)
axis(2, col = "white", at = cumsum(table(annotations$tsneClass)[c(4, 8:11)+1])[-5], labels = NA, col.ticks = "black", tck = 0.04)
axis(4, col = "white", at = cumsum(table(annotations$tsneClass)[c(4, 8:11)+1])[-5], labels = NA, col.ticks = "black", tck = 0.04)
box()


par(mai = c(0.1, 2, 0.1, 1.5))
image(1:13618, 1:1842,log2(ll@expr.data)[,GetClustSample(annotations, cluster = c(1, 0, 3))], 
      col = rev(brewer.pal(11, "RdBu")), 
      breaks = c(-0.884, -0.6, -0.4, -0.25, -0.2, -0.1, 0.1, 0.3, 0.4, 0.45, 0.6, 0.988),
      axes = F, ann = F)
abline(v = cumsum(table(ll@gene_order$chr))[-23], col = "gray60", lwd = 0.8)
abline(h = cumsum(table(annotations$tsneClass)[c(1, 0, 3)+1])[-3], col = "gray60", lwd = 0.8)
axis(1, col = "white", at = cumsum(table(ll@gene_order$chr))[-23], labels = NA, col.ticks = "black", tck = 0.03)
axis(3, col = "white", at = cumsum(table(ll@gene_order$chr))[-23], labels = NA, col.ticks = "black", tck = 0.03)
axis(2, col = "white", at = cumsum(table(annotations$tsneClass)[c(1, 0, 3)+1])[-3], labels = NA, col.ticks = "black", tck = 0.03)
axis(4, col = "white", at = cumsum(table(annotations$tsneClass)[c(1, 0, 3)+1])[-3], labels = NA, col.ticks = "black", tck = 0.03)
box()

par(mai = c(1, 2, 0.1, 1.5))
image(1:13618, 1:1103, log2(ll@expr.data)[,GetClustSample(annotations, cluster = c(2, 5:7))], 
      col = rev(brewer.pal(11, "RdBu")), 
      breaks = c(-0.884, -0.6, -0.4, -0.25, -0.2, -0.1, 0.1, 0.3, 0.4, 0.45, 0.6, 0.988),
      axes = F, ann = F)
abline(v = cumsum(table(ll@gene_order$chr))[-23], col = "gray60", lwd = 0.8)
abline(h = cumsum(table(annotations$tsneClass)[c(2, 5:7)+1])[-4], col = "gray60", lwd = 0.8)
axis(1, col = "white", at = cumsum(table(ll@gene_order$chr))[-23], labels = NA, col.ticks = "black", tck = 0.03)
axis(3, col = "white", at = cumsum(table(ll@gene_order$chr))[-23], labels = NA, col.ticks = "black", tck = 0.03)
axis(2, col = "white", at = cumsum(table(annotations$tsneClass)[c(2, 5:7)+1])[-4], labels = NA, col.ticks = "black", tck = 0.03)
axis(4, col = "white", at = cumsum(table(annotations$tsneClass)[c(2, 5:7)+1])[-4], labels = NA, col.ticks = "black", tck = 0.03)
box()
color.legend(13818, 0, 14288, 860, gradient = "y", rect.col = rev(brewer.pal(11, "RdBu")))
text(c(14438, 14438), c(35, 430, 820), c(-1, 0, 1), srt = 90, xpd = T)
text(15070, 430, "Inferred CNA\n(log2-ration)", srt = 90, xpd = T, adj = 0.5)
y <- c(695.0,1846.5,2701.5,3353.0,3927.5,4594.5,5267.5,5829.5,6332.5,6862.5,7514.0,8270.5,8759.0,9105.0,
       9540.0,10050.0,10749.0,11248.5,11831.5,12495.5,12743.0,12967.5,13372.0)
text(y[-c(13, 15, 18, 21, 22)], rep(-150, 18),c(as.character(c(1:12, 14, 16, 17, 19, 20)), "X"), xpd = T)
par(opar)
dev.off()


annotations_data <- annotations %>% 
  dplyr::mutate(tSNE1 = gliomass[["tsne"]]@cell.embeddings[rownames(annotations), "tSNE_1"],
                tSNE2 = gliomass[["tsne"]]@cell.embeddings[rownames(annotations), "tSNE_2"],
                sample = rownames(annotations)) %>%
  dplyr::arrange(., tsneClass)
P <- ggplot(annotations_data, aes(x = tSNE1, y = tSNE2, colour = tsneClass))


top20 <- gliomass.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC) %>% ungroup(cluster)
featuresgene20 <- c(top20[top20$cluster == "1",]$gene, top20[top20$cluster == "0",]$gene,
                    top20[top20$cluster == "3",]$gene)
top15 <- gliomass.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) %>% ungroup(cluster)
featuresgene15 <- c(top15[top15$cluster == "1",]$gene, top15[top15$cluster == "0",]$gene,
                    top15[top15$cluster == "3",]$gene)
class10_3 <- FindMarkers(gliomass, ident.1 = c(1,0), ident.2 = 3,
                         slot = "data")
class1_03 <- FindMarkers(gliomass, ident.1 = 1, ident.2 = c(0, 3),
                         slot = "data")
class1_sample <- WhichCells(gliomass, idents = 1)
class0_sample <- WhichCells(gliomass, idents = 0)
class3_sample <- WhichCells(gliomass, idents = 3)
class103sample <- gliomass[["RNA"]]@scale.data[, c(class1_sample, class0_sample, class3_sample)]
maker20data <- class103sample[featuresgene20,]

maker20datalong <- maker20data %>% 
  as.data.frame %>%
  rownames_to_column("gene") %>%
  gather(-gene, key = "sample", value = "expression")
maker20datalong$gene <- factor(maker20datalong$gene, levels = featuresgene20)
maker20datalong$sample <- factor(maker20datalong$sample, levels = c(class1_sample, class0_sample, class3_sample))
maker20datalong$expression <- ifelse(maker20datalong$expression > 2.5, 2.5, 
                                     ifelse(maker20datalong$expression < -2.5, -2.5, maker20datalong$expression))
maker20datalong %>%
  ggplot(., aes(x = sample, y = desc(gene))) +
  geom_tile(aes(fill = expression)) +
  scale_fill_gradientn(colors = rev(brewer.pal(11, "RdBu"))) +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())
class10_3_diff <- class10_3 %>% rownames_to_column("gene") %>% filter(abs(avg_logFC)>=1, p_val_adj<0.05)
class1_03_diff <- class1_03 %>% rownames_to_column("gene") %>% filter(abs(avg_logFC)>=1, p_val_adj<0.05)
diffgene <- c(class10_3_diff[class10_3_diff$avg_logFC > 0,]$gene, 
              class1_03_diff[class1_03_diff$avg_logFC<0,]$gene)
diffgenes <- setdiff(diffgene, featuresgene20)

diff <- class103sample[diffgenes,]
diffs <- diff %>%
  as.data.frame %>%
  rownames_to_column("gene") %>%
  gather(-gene, key = "sample", value = "expression")
diffs$gene <- factor(diffs$gene, levels = diffgenes)
diffs$sample <- factor(diffs$sample, levels = c(class1_sample, class0_sample, class3_sample))
diffs$expression <- ifelse(diffs$expression > 2.5, 2.5,
                           ifelse(diffs$expression < -2.5, -2.5, diffs$expression))
diffs %>%
  ggplot(., aes(x = sample, y = gene)) +
  geom_tile(aes(fill = expression)) +
  scale_fill_gradientn(colors = rev(brewer.pal(11, "RdBu"))) +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank())

allgene <- c(featuresgene15, diffgenes15)
alldiff <- class103sample[allgene,]
alldiffs <- alldiff %>%
  as.data.frame %>%
  rownames_to_column("gene") %>%
  gather(-gene, key = "sample", value = "expression")
alldiffs$gene <- factor(alldiffs$gene, levels = rev(allgene))
alldiffs$sample <- factor(alldiffs$sample, levels = c(class1_sample, class0_sample, class3_sample))
alldiffs$expression <- ifelse(alldiffs$expression > 2.5, 2.5,
                           ifelse(alldiffs$expression < -2.5, -2.5, alldiffs$expression))
alldiffs %>%
  ggplot(., aes(x = sample, y = gene)) +
  geom_tile(aes(fill = expression)) +
  scale_fill_gradientn(colors = rev(brewer.pal(11, "RdBu"))) +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank())

diffgenes15 <- setdiff(diffgene, featuresgene15)


ts10 <- FindMarkers(gliomass, ident.1 = 1, ident.2 = 0, min.pct = 0.5, min.cells.feature = 10, logfc.threshold = 0.5) %>%
  rownames_to_column("gene")
ts13 <- FindMarkers(gliomass, ident.1 = 1, ident.2 = 3, min.pct = 0.5, min.cells.feature = 10, logfc.threshold = 0.5) %>%
  rownames_to_column("gene")
ts30 <- FindMarkers(gliomass, ident.1 = 3, ident.2 = 0, min.pct = 0.5, min.cells.feature = 10, logfc.threshold = 0.5) %>%
  rownames_to_column("gene")


ts10_3 <- intersect(ts13[ts13$avg_logFC > 0, ]$gene, ts30[ts30$avg_logFC < 0,]$gene)
ts30_1 <- intersect(ts10[ts10$avg_logFC < 0, ]$gene, ts13[ts13$avg_logFC < 0,]$gene)
ts <- class103sample[c(ts10_3,ts30_1),] %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  gather(-gene, key = "sample", value = "expression") %>%
  mutate(gene = factor(gene, levels = c(ts10_3, ts30_1)),
         sample = factor(sample, levels = c(class1_sample, class0_sample, class3_sample)))

ts$expression <- ifelse(ts$expression > 2, 2,
                         ifelse(ts$expression < -2, -2, ts$expression))
ts %>%
  ggplot(., aes(x = sample, y = gene)) +
  geom_tile(aes(fill = expression)) +
  scale_fill_gradientn(colors = rev(brewer.pal(11, "RdBu"))) +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank())
ts10_3

# ----------------- infer cnv ----------------------------------
refer <- lapply(as.character(c(0, 1, 3, 4, 8, 9, 10, 11)), function(x) {
  ann <- subset(annotations, tsneClass == x)
  n <- nrow(ann)
  set.seed(110)
  index <- sample(x = 1:n, size = as.integer(n / 10), replace = FALSE)
  ann[index,]
})
