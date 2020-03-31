Getgene <- function(gtf) {
  gtfs = read_tsv(gtf, skip = 5, 
                  col_names = FALSE, 
                  col_types = list(X1 = col_character()),
                  n_max = Inf)
  chr = c(1:22, "X", "Y", "MT")
  gtfs = gtfs[gtfs$X1 %in% chr,]
  gtfs = gtfs[gtfs$X3 == "gene",c(1:5, 7, 9)]
  
  extr_gtf = dplyr::mutate(gtfs,
                           ENSEMID = stringr::str_extract(gtfs$X9, pattern = "ENSG\\d{11}"),
                           GENENAME = stringr::str_split(gtfs$X9, ";", simplify = TRUE)[,2] %>% 
                             stringr::str_split(pattern = '\\"', simplify = TRUE) %>% 
                             .[,2])
  
  extr_gtf = dplyr::mutate(extr_gtf,
                           GENENAMES = stringr::str_split(extr_gtf$GENENAME, pattern = "\\.",
                                                          simplify = TRUE)[,1],
                           EDITION = stringr::str_split(extr_gtf$GENENAME, pattern = "\\.",
                                                        simplify = TRUE)[,2])
  extr_gtf$EDITION = ifelse(extr_gtf$EDITION == "", "1", extr_gtf$EDITION)
  extr_gtf$EDITION = as.numeric(extr_gtf$EDITION)
  
  t <- sapply(unique(extr_gtf$GENENAMES), function(x) {da = extr_gtf[extr_gtf$GENENAMES == x,]; da[which.max(da$EDITION),]})
  extr_gtfs <- t(t)
  
  extr_gtf <- apply(extr_gtf, 2, unlist)
  extr_gtfs <- apply(extr_gtfs, 2, unlist)
  
  extr_gtf %<>% data.frame(stringsAsFactors = FALSE) %>%
    dplyr::mutate(X4 = as.numeric(X4),
                  X5 = as.numeric(X5))
  
  extr_gtfs %<>% data.frame(stringsAsFactors = FALSE) %>%
    dplyr::mutate(X4 = as.numeric(X4),
                  X5 = as.numeric(X5))

  result <- list(allgene = extr_gtf,
                 uniquegene = extr_gtfs)
}


GeteCnv <- function(seurat_object,
                    annotations,
                    gtfs,
                    chr = c("X", "Y", "MT"),
                    refgroup,
                    outputdir,
                    scale.data = FALSE) {
  
  chr_index <- c(as.character(1:22), "X", "Y", "MT")
  sort <- data.frame()
  for(i in chr_index) {
    da <- gtfs[gtfs$X1 == i,]
    da <- da[order(da$X4),]
    sort <- rbind.data.frame(sort, da)
  }
  
  infercnv_obj = infercnv::CreateInfercnvObject(raw_counts_matrix=as.data.frame(seurat_object[["RNA"]]@counts),
                                                annotations_file=data.frame(type=paste0("class", annotations$tsneClass),
                                                                            row.names = rownames(annotations)),
                                                delim="\t",
                                                chr_exclude=chr,
                                                gene_order_file=data.frame(chr = sort$X1,
                                                                           start = sort$X4,
                                                                           end = sort$X5,
                                                                           row.names = sort$GENENAMES),
                                                ref_group_names=refgroup)
  
  
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                               out_dir=outputdir, 
                               cluster_by_groups=TRUE,
                               num_threads = 10,
                               png_res = 300,
                               scale_data = scale.data,
                               denoise=TRUE,
                               HMM=FALSE)
}

GetClustSample <- function(annotations = NULL, cluster) {
  if(is.null(annotations)) {
    stop("annotations is empty")
  }
  if( length(cluster) == 0 ){
    stop("cluster is empty")
  }
  result <- data.frame()
  for(i in cluster){
    result <- rbind.data.frame(result, annotations[annotations$tsneClass == i,])
  }
  return(rownames(result))
}
