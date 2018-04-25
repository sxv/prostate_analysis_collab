Sys.setenv(R_MAX_NUM_DLL=5000)

check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

check.biocLites <- function(pkg){
  source("http://bioconductor.org/biocLite.R")
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    biocLite(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

if("slam" %in% rownames(installed.packages()) == FALSE) {
install.packages('devtools')
library(devtools)
slam_url <- "https://cran.r-project.org/src/contrib/Archive/slam/slam_0.1-37.tar.gz"
install_url(slam_url)}

packages<-c("dplyr","tweedie", "statmod", "tightClust", "mclust", "NbClust", "VGAM","ClassDiscovery",  "fitdistrplus", "circlize", "ggbeeswarm","ggplot2", "parallelDist","pheatmap", "magrittr", "MASS", "psych", "RColorBrewer")
biopackages <- c("BiocParallel", "PoiClaClu", "ConsensusClusterPlus", "GenomicFeatures", "DESeq2", "vsn", "sva", "ReportingTools", "IHW","Glimma")
check.packages(packages)
check.biocLites(biopackages)

library(parallelDist)
#library(cluster)
library(circlize)
#library(preprocessCore)
#library(som)
library(dplyr)
library(ConsensusClusterPlus)
library(ggplot2)
library(PoiClaClu)
library("pheatmap")
library("RColorBrewer")
#library("GenomicFeatures")
library("DESeq2")
#library(psych)
library(sva)
library("magrittr")
library(vsn)
library("BiocParallel")
#library("genefilter")
library("ggbeeswarm")
library("IHW")
library("Glimma")
library("ClassDiscovery")
register(MulticoreParam(6))
##----------Functions ------------####
#Generate pheatmap
save_pheatmap_pdf <- function(x, filename, width=50, height=50) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
#Outlier replacement
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.01, .99), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <-NA
  y[x > (qnt[2] + H)] <-NA
  print("The number of genes affected by outliers.")
  na_count <- apply(y, 1, function(x) sum(is.na(x)))
  b <- na_count>0
  print(paste(sum(b), "(", mean(b)*100, "%)"))
  print("Number of genes has max. outliers:")
  print(sum(na_count>=max(na_count)))
  print("max number of outlier per genes:")
  print(max(na_count))
  print("the mean / median of outlier per genes:")
  print(paste0(mean(na_count), " / ", median(na_count)))
  print("The number of genes")
  print("The number of value cells affected by the outlier replacement")
  c_c = is.na(y)
  tb = table(c_c)
  a = tb["FALSE"]
  b = tb["TRUE"]
  value_affected_by_outlier = (b / (a+b)) *100
  value_affected_by_outlier = format(round(value_affected_by_outlier, 2), nsmall = 2)
  print(paste0(value_affected_by_outlier, " %"))
  y <- ifelse(is.na(y), rowMeans(x), y)
  y
}
#Geometric mean
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
#Get the result with lowest changes
get_recursive_cluster_number <- function(x=NA) {
  #x = x[order(x$k),]
  length_x = nrow(x)
  select_to_cluster = -1
  course_changed= F
  #print(paste0("Mean: ", mean(x$mean_vle)))
  #print(paste0("Median:", median(x$mean_vle)))
  #print(gm_mean(x$mean_vle))
  d_list = as.matrix(dist(x$mean_vle))[,1]
  #print(mean(d_list))
  #print(median(d_list))
  #print(gm_mean(d_list))
  #print(d_list)
  x_select = x
  x_select$b = c(rev(d_list))
  #print(x_select)
  x_select = x_select[x_select$b>mean(d_list),]
  #print(x_select)
  x_select = x_select[order(x_select$b),]
  #print(x_select)
  select_to_cluster = x_select$k[1]
  s_index = match(select_to_cluster,x$k)
  s_index = s_index[1]
  if (s_index <length_x & s_index>1) {
    previous_mean_itm = x$mean_vle[s_index+1]
    previous_sem_itm = x$sem_vle[s_index+1]
    current_mean_itm = x$mean_vle[s_index]
    current_sem_itm = x$sem_vle[s_index]
    next_mean_itm = x$mean_vle[s_index-1]
    next_sem_itm = x$sem_vle[s_index-1]
    print(select_to_cluster)
    if ((previous_mean_itm-previous_sem_itm)<=(current_mean_itm+current_sem_itm)) {
      select_to_cluster = x$k[s_index+1]
      print("Changed back")
    }
    if ((current_mean_itm-current_sem_itm) < (next_mean_itm+next_sem_itm)) {
      if ((next_mean_itm+next_sem_itm)<(previous_mean_itm-previous_sem_itm)) {
        select_to_cluster = x$k[s_index-1]
        print("Changed forward")
      }
    }
  }
  return(select_to_cluster)
}
#Get the possilbe cluster number from COCA
get_cluster_number <- function(icl=NA, prefix_id, normalization_id, data_type_id, k_max=20) {
  #Identify the best cluster for samples base on item cluster (lower value is better)
  results_icl_mean = list()
  results_icl_sem = list()
  results_icl_k = list()
  results_sig_chnages = list()
  select_k_size_sem_selected = 0
  previous_k_mean = 0
  for (counter in c(1:k_max)) {
    list_of_vl = icl$itemConsensus[icl$itemConsensus$k==counter,]
    list_of_vl =na.omit(list_of_vl)
    m_v = mean(list_of_vl$itemConsensus)
    sem_v = sd(list_of_vl$itemConsensus)/sqrt(length(list_of_vl$itemConsensus))
    if (counter<k_max & !is.na(m_v)) {
      list_of_vl = icl$itemConsensus[icl$itemConsensus$k==counter+1,]
      list_of_vl =na.omit(list_of_vl)
      m_v_2 = mean(list_of_vl$itemConsensus)
      sem_v_2 = sd(list_of_vl$itemConsensus)/sqrt(length(list_of_vl$itemConsensus))
      if ((m_v-sem_v) >= (m_v_2+sem_v_2)) {
        results_icl_k[counter+1] = counter+1
        results_icl_mean[counter+1] = m_v_2
        results_icl_sem[counter+1] = sem_v_2
      }
    }
  }

  table_icl = data.frame(k=unlist(results_icl_k), mean_vle = unlist(results_icl_mean), sem_vle = unlist(results_icl_sem))
  table_icl = table_icl[order(table_icl$mean_vle),]
  print(paste0(prefix_id, "_", normalization_id, "_", data_type_id, ":"))
  print(table_icl)

  recommend_cluster = get_recursive_cluster_number(table_icl)#table_icl$k[1]
  print(paste("Cluster recommended:", recommend_cluster))
  png(paste0(prefix_id, "_", normalization_id, "_", data_type_id,"_cluster_sample.png"))
  plot(x=table_icl$k, y=table_icl$mean_vle)
  dev.off()
  return(recommend_cluster)
}
#Get the possible cluster number based on consensus
get_cluster_number_based_on_Consensus <- function(icl=NA, prefix_id, normalization_id, data_type_id) {
  tabl = data.frame(icl$clusterConsensus)
  table_m = list()
  table_k = list()
  counter = 1
  tabl = na.omit(tabl)
  for (x in c(1:max(tabl$k))) {
    table_m[counter] = mean(tabl$clusterConsensus[tabl$k==x])
    table_k[counter] = x
  }
  x_result <- tabl %>% group_by(k) %>%
    summarise(clusterConsensus = mean(clusterConsensus)) %>% arrange(desc(clusterConsensus))# %>% filter(clusterConsensus<=ic95$conf.int[2] & clusterConsensus>=ic95$conf.int[1])
  x_result = na.omit(x_result)

  print(paste0(prefix_id, "_", normalization_id, "_", data_type_id, ":"))
  print(x_result)
  k_sample = 1
  if (nrow(x_result)>0) {
    counter = 1

    previous_value = 0
    for (x_k in signif(x_result$clusterConsensus,3)) {
      if (x_k>=previous_value) {
        previous_value = x_k
        k_sample = x_result$k[counter]
      } else {
        break
      }
      counter = counter +1
    }
  }
  return(k_sample)
}
#Get the tight cluster
getTightCluster <- function(x=NA, stop_cut_off=0.5, target_max=20, standardize.gene=F, Start_with_max=T) {
  incremental_increase = 1

  target_k = 2
  k_min = target_k +5
  tclust2 = NA
  if (Start_with_max){

    tclust2<-tight.clust(x, target=target_max, k.min=target_max+5, random.seed=12345, standardize.gene=standardize.gene)
    print(table(tclust2$cluster))
  } else{
  repeat{
    tclust2<-tight.clust(x, target=target_k, k.min=k_min, random.seed=12345, standardize.gene=standardize.gene)
    c_r = length(tclust2$cluster[tclust2$cluster==-1]) / length(tclust2$cluster)
    target_k = target_k + incremental_increase
    k_min = target_k + 5
    if(c_r<=stop_cut_off | target_k >=target_max){
      print(table(tclust2$cluster))
      break
    }
  }
  }
  return(tclust2)
}
#Get the heatmap from tight cluster
#Order sample either by hclust or the columne name in col_anno: SampleType
generate_pheatmap_from_tight_cluster <- function (x, standardize.gene = T, col_anno=NULL, order.sample = "hclust", plot.noise=F) {
  data <- x$data
  cluster <- x$cluster
  size <- x$size
  data_t = data[!duplicated(rownames(data),fromLast = T),]
  cluster = cluster[!duplicated(rownames(data),fromLast = T)]
  data = data_t

  if (standardize.gene)
    data <- t(scale(t(data)))
  l <- 1:max(cluster)
  if (plot.noise)
    l <- c(l, -1)
  order.genes <- unlist(lapply(l, function(i) which(cluster == i)))
  data.plot <- data[order.genes, ]
  cluster_ordered <- cluster[order.genes]
  SampleType = NA
  if (order.sample=="hclust") {
    order.sample = T
  } else {
    l_x <- levels(col_anno[,order.sample])
    sample_list <- col_anno[,order.sample]
    order.samples <- unlist(lapply(l_x, function(i) which(sample_list == i)))
    #print(length(order.samples))
    #print(dim(data.plot))
    SampleType = order.sample
    data.plot <- data.plot[, order.samples]
    #print(dim(data.plot))
    #print("Done")s
    order.sample = F
  }

  col = colorRampPalette(c("green", "black", "red"))(21)
  info_row = data.frame('gene cluster' = as.factor(cluster_ordered[cluster_ordered %in% l]))
  #print(dim(info_row))
  rownames(info_row) = rownames(data.plot)
  ann_colors = list('gene cluster' = colorRampPalette(c("red", "black"))(max(cluster)))
  names(ann_colors$'gene cluster') = levels(info_row$gene.cluster)
  if (!is.null(col_anno)) {
    ann_colors$SampleType = rainbow(length(unique(col_anno[,SampleType])))
    names(ann_colors$SampleType) = unique(col_anno$SampleType)
  }

  pt <- pheatmap(data.plot, annotation_row = info_row, annotation_col = col_anno,  color = col, clustering_distance_cols="correlation",
                 annotation_colors = ann_colors,
                 cluster_cols=order.sample,
                 cluster_rows = F,
                 fontsize = 6.5,
                 fontsize_row=1,
                 fontsize_col = 2)
  return(pt)
}

##--------DESEQ-------------------######
#Unmark these area for debuging. Including required info for Debug (Pls run the part ofRNA_Analysis.R before running this command)
#x=gene_data_raw_for_DESEQ2;
#initial_filter = initial_filter
#set_cutoff = cutoff_level
#Col_Data=sample_Data_for_seq
#prefix=factor_X
#Run_Id=counter_id
#description_txt=factor_X
#PlayWithResults=T
#transform_regulator="vst"
#group_factor=c("sample_type", "batch", "SampleID")
#contrast_list=c(levels(sample_Data_for_seq$sample_type))
#Qutlier_type= None, Quantile,  MVO
RunSEQ <- function(x=NULL, dds = NULL, initial_filter="rowSum", Col_Data=NULL, outlier_type="Quantile", transform_regulator="vst", set_cutoff=1, Run_Id=0, description_txt="", prefix="X", group_factor=c("sample_type", "batch", "SampleID"),  contrast_list=c(), PlayWithResults=FALSE) {
  filn_n = paste0(prefix, "_results.RData")
  if (file.exists(filn_n)) {
  load(file =  filn_n)
    }
  print(paste0("Run id: ",Run_Id))
  print("++++++++++++++++++++++++")
  print(description_txt)
  print("++++++++++++++++++++++++")
  print("Proc: Loading the data...")
  if (!is.null(x)) {
    dds <- DESeqDataSetFromMatrix(countData = x, colData = Col_Data, design= ~ sample_type+batch)
  } else {
    if (is.null(dds)) { stop("No values defined please define all parameter and use either x or dds. x will be preferred.")}
  }

  print("Done: Loading the data...")
  Save_row =nrow(dds)
  if (initial_filter=="rowSum") {
    print("Runing the filter...")
    dds <- dds[ rowSums(counts(dds)) > set_cutoff, ]
    print(paste0("gene number after the row sum filter: ", nrow(dds), " reduced by ", Save_row-nrow(dds), " from ", Save_row))
  } else {
    print("No filter...")
  }
  df <- as.data.frame(colData(dds)[,group_factor])
  notAllZero <- (rowSums(counts(dds))>0)
  if (transform_regulator=="rlog" | transform_regulator=="all") {
    print("Runing transformation based on rlog...")
    ##----rlog----------------------------------------------------------------
    rld <- rlog(dds, blind = FALSE)
    if (transform_regulator=="rlog") {
      vsd <- rld
    }
    msd <- meanSdPlot(assay(rld[notAllZero,]))
    msd$gg
    ggsave(paste0(prefix,"_plot_rlt_effect_of_transformation_on_the_variance.png"), width = 5, height = 5, unit="in", limitsize=F)

    ##----rlog-sample-heatmap----------------------------------------------------------
    sampleDists <- dist(t(assay(rld)))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- paste(df$SampleID, df$batch, sep="-")
    colnames(sampleDistMatrix) <- NULL
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    ptl_heatmap <- pheatmap(sampleDistMatrix,
             clustering_distance_rows=sampleDists,
             clustering_distance_cols=sampleDists,
             col=colors)
    save_pheatmap_pdf(ptl_heatmap, paste0(prefix,"_sample_rld_heatmap.pdf"))
  }
  save.image(file =  paste0(prefix, "_results.RData"))
  if (transform_regulator=="vst" | transform_regulator=="none" | transform_regulator=="all") {
    ## ----vst-----------------------------------------------------------------
    print("Runing transformation based on vst...")
    vsd <- vst(dds, blind = FALSE)
    if (transform_regulator == "vst" | transform_regulator=="none")
    {
      rld <- vsd
    }
    med <- meanSdPlot(assay(vsd[notAllZero,]))
    med$gg
    ggsave(paste0(prefix,"_plot_vst_effect_of_transformation_on_the_variance.png"), width = 5, height = 5, unit="in", limitsize=F)

    sampleDists <- dist(t(assay(vsd)))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- paste(df$SampleID, df$batch, sep="-")
    colnames(sampleDistMatrix) <- NULL
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    ptl_heatmap <- pheatmap(sampleDistMatrix,
             clustering_distance_rows=sampleDists,
             clustering_distance_cols=sampleDists,
             col=colors)
    save_pheatmap_pdf(ptl_heatmap, paste0(prefix,"_sample_vst_heatmap.pdf"))
  }
  save.image(file =  paste0(prefix, "_results.RData"))
  head(assay(rld), 3)
  head(assay(vsd), 3)
  ## ------------------------------------------------------------------------
  print("Generate the distance martix for samples.")
  sampleDists <- dist(t(assay(rld)))
  #sampleDists

  ## ----distheatmap, fig.width = 6.1, fig.height = 4.5----------------------
  print("Generate a heat map based on distance martix for samples...")
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(rld$SampleID , rld$batch, sep = " - " )
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  ptl_heatmap = pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors)
  save_pheatmap_pdf(ptl_heatmap, paste0(prefix,"_distance_heatmap.pdf"))
  save.image(file =  paste0(prefix, "_results.RData"))
  ## ------------------------------------------------------------------------
  library("PoiClaClu")
  print("Generate a heat map based on poisson distance martix for samples...")
  poisd <- PoissonDistance(t(counts(dds)))

  ## ----poisdistheatmap, fig.width = 6.1, fig.height = 4.5------------------
  samplePoisDistMatrix <- as.matrix( poisd$dd )
  rownames(samplePoisDistMatrix) <- paste(Col_Data$SampleID , rld$batch, sep=" - " )
  colnames(samplePoisDistMatrix) <- NULL
  ptl_heatmap = pheatmap(samplePoisDistMatrix,
           clustering_distance_rows = poisd$dd,
           clustering_distance_cols = poisd$dd,
           col = colors)
  save_pheatmap_pdf(ptl_heatmap, paste0(prefix,"_poisson_distance_heatmap.pdf"))
  save.image(file =  paste0(prefix, "_results.RData"))
  ## ----plotpca, fig.width=6, fig.height=4.5--------------------------------
  print("Generate a PCA plot...")
  png(paste0(prefix,"_PCA_sample_type_batch.png"))
  plotPCA(rld, intgroup = c("sample_type", "batch"))
  dev.off()
  ## ------------------------------------------------------------------------
  print("Generate PCA plot...")
  pcaData <- plotPCA(rld, intgroup = c( "SampleID", "batch"), returnData = TRUE)
  #pcaData
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ## ----ggplotpca, fig.width=6, fig.height=4.5------------------------------
  ggplot(pcaData, aes(x = PC1, y = PC2, color = SampleID, shape = batch)) +
    geom_point(size =3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed()
  ggsave(paste0(prefix,"_PCA1_2_sample_type_batch.png"), width = 50, height = 50, unit="in", limitsize=F)
  ## ----mdsrlog, fig.width=6, fig.height=4.5--------------------------------
  print("Generate MDS plot...")
  mds <- as.data.frame(colData(rld))  %>%
    cbind(cmdscale(sampleDistMatrix))
  ggplot(mds, aes(x = `1`, y = `2`, color = SampleID, shape = batch)) +
    geom_point(size = 3) + coord_fixed()
  ggsave(paste0(prefix,"_mds_sample_type_batch.png"), width = 50, height = 50, unit="in", limitsize=F)
  ## ----mdspois, fig.width=6, fig.height=4.5--------------------------------
  print("Generate poisson MDS plot...")
  mdsPois <- as.data.frame(colData(dds)) %>%
    cbind(cmdscale(samplePoisDistMatrix))
  ggplot(mdsPois, aes(x = `1`, y = `2`, color = SampleID, shape = batch)) +
    geom_point(size = 3) + coord_fixed()
  ggsave(paste0(prefix,"_mds_poisson_sample_type_batch.png"), width = 50, height = 50, unit="in", limitsize=F)
  ## ----RunDESEQ------------------------------------------------------------
  Y = list()
  Y$dds <- DESeq(dds, parallel = T)
  Y$res <- results(Y$dds, parallel = T)
  save.image(file =  paste0(prefix, "_results.RData"))
  ## ---RunDESEQ_with_LRT----------------------------------------------------
  dds_lrt <- DESeq(dds, test="LRT", reduced=~batch, parallel = T)
  res_lrt <- results(dds_lrt, parallel = T)
  Y$dds_lrt = dds_lrt
  Y$res_lrt = res_lrt
  save.image(file =  paste0(prefix, "_results.RData"))
  ##--------------------------------------------------------------------------------
  print("Generate transformed counts")
  print(paste0(transform_regulator, " approach is selected."))
  ##prepare the parameters for the next step...
  notAllZero <- (rowSums(counts(Y$dds))>0)
  select <- order(rowMeans(counts(Y$dds,normalized=TRUE)),
                  decreasing=TRUE)
  nt <- normTransform(Y$dds) # defaults to log2(x+1)
  log2.norm.counts <- assay(nt)[select,]
  ##----rlog/vsd-gene-heatmap------------------------------------------------------------
  ptl_heatmap = pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
           cluster_cols=FALSE, annotation_col=df)
  save_pheatmap_pdf(ptl_heatmap, paste0(prefix,"_gene_rld_heatmap.pdf"))
  ## ----
  ptl_heatmap = pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
                         cluster_cols=FALSE, annotation_col=df)
  save_pheatmap_pdf(ptl_heatmap, paste0(prefix,"_gene_vst_2_heatmap.pdf"))
  ## ----rldplot, fig.width = 6, fig.height = 2.5----------------------------
  print("Generate Scatterplot of transformed counts from first two samples.")
  dds_f <- Y$dds
  dds_f <- estimateSizeFactors(dds)
  ###------------
  df_x <- bind_rows(
    as_data_frame(log2(counts(dds_f, normalized=TRUE)+1)) %>%
      mutate(transformation = "log2(x + 1)"),
    as_data_frame(assay(rld)) %>% mutate(transformation = transform_regulator))

  colnames(df_x)[1:2] <- c("x", "y")
  ggplot(df_x, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)
  ggsave(paste0(prefix,"_Scatter_plot_transformd_count_from_samples.png"))

  ## ----Count outlier detection from DES----------------------------
  W <- Y$res$stat
  maxCooks <- apply(assays(Y$dds)[["cooks"]],1,max)
  idx <- !is.na(W)
  png(paste0(prefix,"_Count_outlier_detection_from_DES.png"))
  plot(rank(W[idx]), maxCooks[idx], xlab="rank of Wald statistic",
       ylab="maximum Cook's distance per gene", main="DES",
       ylim=c(0,5), cex=.4, col=rgb(0,0,0,.3))
  m <- ncol(Y$dds)
  p <- 3
  abline(h=qf(.99, p, m - p))
  dev.off()
  ## ----Count outlier detection from LRT ----------------------------
  W <- Y$res_lrt$stat
  maxCooks <- apply(assays(Y$dds_lrt)[["cooks"]],1,max)
  idx <- !is.na(W)
  png(paste0(prefix,"_Count_outlier_detection_from_LRT.png"))
  plot(rank(W[idx]), maxCooks[idx], xlab="rank of Wald statistic",
       ylab="maximum Cook's distance per gene", main="LRT",
       ylim=c(0,5), cex=.4, col=rgb(0,0,0,.3))
  m <- ncol(Y$dds_lrt)
  p <- 3
  abline(h=qf(.99, p, m - p))
  dev.off()
  ## --- Filtering criteria LRT----
  png(paste0(prefix,"_FilteringCriteria_LRT.png"))
  plot(Y$res_lrt$baseMean+1, -log10(Y$res_lrt$pvalue),
       log="x", xlab="mean of normalized counts", main="DES",
       ylab=expression(-log[10](pvalue)),
       ylim=c(0,30),
       cex=.4, col=rgb(0,0,0,.3))
  dev.off()
  ## --- Filtering criteria DES----
  png(paste0(prefix,"_FilteringCriteria_DES.png"))
  plot(Y$res$baseMean+1, -log10(Y$res$pvalue),
       log="x", xlab="mean of normalized counts", main="DES",
       ylab=expression(-log[10](pvalue)),
       ylim=c(0,30),
       cex=.4, col=rgb(0,0,0,.3))
  dev.off()
  ## --- Save RData for later use ----
  print("Save the memory image to the hard disk...")
  save.image(file =  paste0(prefix, "_results.RData"))

  ## --- how the filtering ameliorates the multiple testing problem
  png(paste0(prefix,"_Histogram_FilteringIssue_Pvalue_cutoff_DES.png"))
  use <- Y$res$baseMean > metadata(Y$res)$filterThreshold
  h1 <- hist(Y$res$pvalue[!use], breaks=0:50/50, plot=FALSE)
  h2 <- hist(Y$res$pvalue[use], breaks=0:50/50, plot=FALSE)
  colori <- c(`do not pass`="khaki", `pass`="powderblue")

  barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
          col = colori, space = 0, main = "DES", ylab="frequency")
  text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
       adj = c(0.5,1.7), xpd=NA)
  legend("topright", fill=rev(colori), legend=rev(names(colori)))
  dev.off()
  ## --- how the filtering ameliorates the multiple testing problem
  png(paste0(prefix,"_Histogram_FilteringIssue_Pvalue_cutoff_LRT.png"))
  use <- Y$res_lrt$baseMean > metadata(Y$res_lrt)$filterThreshold
  h1 <- hist(Y$res_lrt$pvalue[!use], breaks=0:50/50, plot=FALSE)
  h2 <- hist(Y$res_lrt$pvalue[use], breaks=0:50/50, plot=FALSE)
  colori <- c(`do not pass`="khaki", `pass`="powderblue")
  barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
          col = colori, space = 0, main = "LRT", ylab="frequency")
  text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
       adj = c(0.5,1.7), xpd=NA)
  legend("topright", fill=rev(colori), legend=rev(names(colori)))
  dev.off()

  ## ---- Generate a heatmap for gene expression and sample (Perform batch effect reduction)----

  ## ---- DESEQ data (Non-zeros values) ----
  # ---- 1. Normalize ---
  counter_raw_data <- Y$dds[ rowSums(counts(Y$dds)) > 1, ]

  sampledata <- colData(Y$dds)

  x_type = "VST_ZScore"
  for (x_type in c("VST", "VST_ZScore", "median", "VST_RANK_Row", "VST_RANK_Col", "DESEQ", "log2", "RANK_Row","RANK_Col", "RANK_Row_DESEQ","RANK_Col_DESEQ", "log2_DESEQ", "ZScore","ZScore_DESEQ")) {
    ## Selection of normalisation approaches #####
    print(paste0("Selected normalisation approach: ", x_type))
    if (x_type == "VST") {
      counter_dd = assay(vst(counter_raw_data))
    }
    if (x_type == "DESEQ") {
      counter_dd = counts(counter_raw_data,normalized=TRUE)
    }
    if (x_type == "log2_DESEQ") {
      counter_dd = log2(counts(counter_raw_data,normalized=TRUE)+1)
    }
    if (x_type == "log2") {
      counter_dd = log2(counts(counter_raw_data,normalized=FALSE)+1)
    }
    if (x_type== "RANK_Row") {
      counter_dd = counts(counter_raw_data,normalized=FALSE)
      data_seq_row_normalize_quantile = counter_dd
      for (i in 1:nrow(data_seq_row_normalize_quantile)) {
        q= unlist(counter_dd[i,])
        perc.rank <- ecdf(q)
        data_seq_row_normalize_quantile[i,] = perc.rank(counter_dd[i,])
      }
      counter_dd = data_seq_row_normalize_quantile
    }
    if (x_type== "RANK_Col") {
      counter_dd = counts(counter_raw_data,normalized=FALSE)
      data_seq_row_normalize_quantile = counter_dd
      for (i in 1:nrow(data_seq_row_normalize_quantile)) {
        q= unlist(counter_dd[,i])
        perc.rank <- ecdf(q)
        data_seq_row_normalize_quantile[,i] = perc.rank(counter_dd[,i])
      }
      counter_dd = data_seq_row_normalize_quantile
    }
    if (x_type== "RANK_Col_DESEQ") {
      counter_dd = counts(counter_raw_data,normalized=TRUE)
      data_seq_row_normalize_quantile = counter_dd
      for (i in 1:nrow(data_seq_row_normalize_quantile)) {
        q= unlist(counter_dd[,i])
        perc.rank <- ecdf(q)
        data_seq_row_normalize_quantile[,i] = perc.rank(counter_dd[,i])
      }
      counter_dd = data_seq_row_normalize_quantile
    }
    if (x_type== "RANK_Row_DESEQ") {
      counter_dd = counts(counter_raw_data,normalized=TRUE)
      data_seq_row_normalize_quantile = counter_dd
      for (i in 1:nrow(data_seq_row_normalize_quantile)) {
        q= unlist(counter_dd[i,])
        perc.rank <- ecdf(q)
        data_seq_row_normalize_quantile[i,] = perc.rank(counter_dd[i,])
      }
      counter_dd = data_seq_row_normalize_quantile
    }
    if (x_type== "ZScore") {
      counter_dd = counts(counter_raw_data,normalized=F)
      counter_dd = scale(t(counter_dd), center = T, scale = T)
      counter_dd = t(counter_dd)
    }
    if (x_type== "ZScore_DESEQ") {
      counter_dd = counts(counter_raw_data,normalized=T)
      counter_dd = scale(t(counter_dd), center = T, scale = T)
      counter_dd = t(counter_dd)
    }
    if (x_type=="VST_ZScore") {
      counter_dd = assay(vst(counter_raw_data))
      counter_dd = scale(t(counter_dd), center = T, scale = T)
      counter_dd = t(counter_dd)
      table(is.na(counter_dd))
    }
    if (x_type=="VST_RANK_Row") {
      counter_dd = assay(vst(counter_raw_data))
      data_seq_row_normalize_quantile = counter_dd
      for (i in 1:nrow(data_seq_row_normalize_quantile)) {
        q= unlist(counter_dd[i,])
        perc.rank <- ecdf(q)
        data_seq_row_normalize_quantile[i,] = perc.rank(counter_dd[i,])
      }
      counter_dd = data_seq_row_normalize_quantile
    }
    if (x_type=="VST_RANK_Col") {
      counter_dd = assay(vst(counter_raw_data))
      data_seq_row_normalize_quantile = counter_dd
      for (i in 1:nrow(data_seq_row_normalize_quantile)) {
        q= unlist(counter_dd[,i])
        perc.rank <- ecdf(q)
        data_seq_row_normalize_quantile[,i] = perc.rank(counter_dd[,i])
      }
      counter_dd = data_seq_row_normalize_quantile
      #plot(counter_dd)
    }
    if (x_type=="median") {
      data_seq_median = scale(t(counter_dd),scale=apply(counter_dd,1,mad))
      counter_dd = t(data_seq_median)
    }
    #Remove outlier
    ##Show the distribution of the count of raw exp data (Before correction and outlier replacement by mean)####
    png(paste0(prefix,"_", x_type, "_distribution_before_outlier_replacment_by_mean_and_batch_effect_correction_gene_expression_data.png"))
    colramp = colorRampPalette(c(3,"white",2))(ncol(Y$dds))
    ptx = plot(density(counter_dd[,1]),col=colramp[1],lwd=3,ylim=c(0,10), xlim=c(min(counter_dd),max(counter_dd)))
    for(i in 2:ncol(Y$dds)){lines(density(counter_dd[,i]),lwd=3,col=colramp[i])}
    dev.off()
    png(paste0(prefix,"_", x_type, "_m_gm_median_distribution_before_outlier_replacment_by_mean_and_batch_effect_correction_gene_expression_data.png"))
    ptx = plot(density(counter_dd[,1]),col=colramp[1],lwd=3,ylim=c(0,10), xlim=c(min(counter_dd),max(counter_dd)))
    for(i in 2:ncol(Y$dds)){lines(density(counter_dd[,i]),lwd=3,col=colramp[i])}
    lines(density(apply(counter_dd, 1, geometric.mean)), lwd=3, col="black")
    lines(density(rowMeans(counter_dd)), lwd=3, col="blue")
    lines(density(apply(counter_dd, 1, median)), lwd=3, col="Yellow")
    dev.off()

    print("dim before outlier replacement by mean 1%-99% margin")
    print(dim(counter_dd))
    if (outlier_type=="Quantile") {
      counter_dd <- remove_outliers(counter_dd)
      counter_dd <- na.omit(counter_dd)
    }
    else if (outlier_type == "MVO") {
      print("Not implemented.")
    }
    else {
      print("Additional outlier not performed.")
    }

    print("dim after outlier rep")
    print(dim(counter_dd))
    table(is.na(counter_dd))
    ##Show the distribution of the count of raw exp data (After correction and outlier replacement by mean)####
    png(paste0(prefix,"_", x_type, "_distribution_before_outlier_replacment_by_mean_and_batch_effect_correction_gene_expression_data.png"))
    colramp = colorRampPalette(c(3,"white",2))(ncol(Y$dds))
    plot(density(counter_dd[,1]),col=colramp[1],lwd=3,ylim=c(0,10), xlim=c(min(counter_dd),max(counter_dd)))
    for(i in 2:ncol(Y$dds)){lines(density(counter_dd[,i]),lwd=3,col=colramp[i])}
    #lines(density(apply(counter_dd, 1, geometric.mean)), lwd=3, col="black")
    lines(density(rowMeans(counter_dd)), lwd=3, col="blue")
    lines(density(apply(counter_dd, 1, median)), lwd=3, col="Yellow")
    dev.off()
    print(dim(counter_dd))

    ## Remove batch effect ####
    mod = model.matrix(~as.factor(sample_type) + as.factor(batch2), data=sampledata)
    mod0 = model.matrix(~+ as.factor(batch2), data=sampledata)
    modcombat = model.matrix(~as.factor(sample_type), data=sampledata)
    gene_exp_corrected = ComBat(dat=counter_dd, batch=sampledata$batch2, mod=modcombat, par.prior=TRUE)
    print(paste0("The dimension of the corrected gene matrix: ",dim(gene_exp_corrected)[1], "x", dim(gene_exp_corrected)[2]))
    pValuesComBat = f.pvalue(gene_exp_corrected,mod,mod0)
    qValuesComBat = p.adjust(pValuesComBat,method="BH")
    gene_list_of_interest = names(qValuesComBat[qValuesComBat<0.05])
    gene_list_of_interest = gene_exp_corrected[gene_list_of_interest,]
    print("The dimension of varied genes not associated with batch effects...")
    print(dim(gene_list_of_interest))

    ##Show the distribution of the count of corrected exp data ####
    png(paste0(prefix,"_", x_type, "_distribution_of_normalized_gene_expression_data.png"))
    colramp = colorRampPalette(c(3,"white",2))(ncol(Y$dds))
    plot(density(gene_exp_corrected[,1]),col=colramp[1],lwd=3,ylim=c(0,3), xlim=c(min(gene_exp_corrected),max(gene_exp_corrected)))
    for(i in 2:ncol(Y$dds)){lines(density(gene_exp_corrected[,i]),lwd=3,col=colramp[i])}
    lines(density(rowMeans(gene_exp_corrected)), lwd=3, col="blue")
    lines(density(apply(gene_exp_corrected, 1, median)), lwd=3, col="Yellow")
    dev.off()
    ##Show the distribution plot of the count of varied genes ----
    png(paste0(prefix, "_", x_type,"_selected_genes_of_interest_distribution_of_normalized_gene_expression_data.png"))
    colramp = colorRampPalette(c(3,"white",2))(ncol(Y$dds))
    plot(density(gene_list_of_interest[,1]),col=colramp[1],lwd=3,ylim=c(0,3), xlim=c(min(gene_list_of_interest), max(gene_list_of_interest)))
    for(i in 2:ncol(Y$dds)){lines(density(gene_list_of_interest[,i]),lwd=3,col=colramp[i])}
    lines(density(apply(gene_list_of_interest, 1, geometric.mean)), lwd=3, col="black")
    lines(density(rowMeans(gene_list_of_interest)), lwd=3, col="blue")
    lines(density(apply(gene_list_of_interest, 1, median)), lwd=3, col="Yellow")
    dev.off()

    ##Identfiy the distribution ####
    library(fitdistrplus)
    #--------- Corrected gene data ----------------#
    png(paste0(prefix, "_", x_type ,"_distribution_of_batch_corrected_gene_expression_data.png"))
    x_value = as.numeric(gene_exp_corrected)
    print(descdist(x_value))
    print(t.test(x_value))
    dev.off()
    #------- Gene of interests ----------------#
    png(paste0(prefix, "_", x_type ,"_distribution_of_selected_genes.png"))
    x_value = as.numeric(gene_list_of_interest)
    print(descdist(x_value))
    print(t.test(x_value))
    dev.off()
    #------ GET CDF  ------------------------------#
    #--------- Corrected gene data ----------------#
    png(paste0(prefix, "_", x_type ,"_CDF_all_corrected_gene_of_interest.png"))
    P = ecdf(gene_exp_corrected)
    plot(P)
    dev.off()
    #------- Gene of interests ----------------#
    png(paste0(prefix, "_", x_type ,"_CDF_gene_of_interest.png"))
    P = ecdf(gene_list_of_interest)
    plot(P)
    dev.off()

    ##---- Fit to different distribution (Gamma, Beta and lognorm) ----#####
    library(VGAM)
    if (min(gene_list_of_interest)>=0) {
      #-------- Gamma -------#
      print("Fitting data to gamma distribution")
      png(paste0(prefix,"_", x_type ,"_gene_of_interest_fitting_to_gamma_distribution.png"))
      fit.gamma <- fitdist(as.numeric(gene_list_of_interest),distr="gamma")
      plot(fit.gamma)
      dev.off()
      png(paste0(prefix,"_", x_type ,"_corrected_genes_fitting_to_gamma_distribution.png"))
      fit.gamma <- fitdist(as.numeric(gene_exp_corrected),distr="gamma")
      plot(fit.gamma)
      dev.off()
      #-------- Beta distribution ------#
      print("Fitting data to beta distribution")
      png(paste0(prefix, "_", x_type, "_gene_of_interest_fitting_to_Beta_distribution.png"))
      fit.beta <- fitdist(as.numeric(gene_list_of_interest),distr=dbetanorm,start=list(shape1=1,shape2=1), fix.arg=list(mean=0,sd=1)) #fitdist(as.numeric(data_seq_raw_normalize_quantile), "beta", lower = c(0, 0))
      plot(fit.beta)
      dev.off()
      png(paste0(prefix, "_", x_type, "_corrected_genes_fitting_to_Beta_distribution.png"))
      fit.beta <- fitdist(as.numeric(gene_exp_corrected),distr=dbetanorm,start=list(shape1=1,shape2=1), fix.arg=list(mean=0,sd=1)) #fitdist(as.numeric(data_seq_raw_normalize_quantile), "beta", lower = c(0, 0))
      plot(fit.beta)
      dev.off()
      #------ log normal distribution -----#
      print("Fitting data to log-normal distribution")
      png(paste0(prefix,"_", x_type ,"_gene_of_interest_fitting_to_log-normal_distribution.png"))
      fit.lognormal <- params <- fitdist(as.numeric(gene_list_of_interest),distr="lnorm") #fitdist(as.numeric(data_seq_raw_normalize_quantile), "beta", lower = c(0, 0))
      plot(fit.lognormal)
      dev.off()
      png(paste0(prefix,"_", x_type , "_corrected_genes_fitting_to_log-normal_distribution.png"))
      fit.lognormal <- params <- fitdist(as.numeric(gene_exp_corrected),distr="lnorm") #fitdist(as.numeric(data_seq_raw_normalize_quantile), "beta", lower = c(0, 0))
      plot(fit.lognormal)
      dev.off()
    }

    #------ uniform distribution -----#
    print("Fitting data to uniform distribution")
    png(paste0(prefix,"_", x_type ,"_gene_of_interest_fitting_to_unif_distribution.png"))
    fit.unif <- fitdist(as.numeric(gene_list_of_interest),distr="unif")
    plot(fit.unif)
    dev.off()
    png(paste0(prefix,"_", x_type ,"_corrected_genes_fitting_to_unif_distribution.png"))
    fit.unif <- fitdist(as.numeric(gene_exp_corrected),distr="unif")
    plot(fit.unif)
    dev.off()
  ####---- Deactivated --------############
    #----- Tweedie distribution ----#
    #library(tweedie)
    #library(statmod)
    #png(paste0(prefix,"_", x_type ,"_corrected_genes_fitting_to_unif_distribution.png"))
    #dis_tweedie = tweedie.profile(as.numeric(gene_exp_corrected)~1, do.plot=FALSE)
    #print(dx)[c('xi.max', 'phi.max')]
    #dev.off()
    #---- Calculate th poisson distance for samples ----##
    #poisd_sample_gene_of_interest <- PoissonDistance(t(gene_list_of_interest))$dd
    #Check if this has NA
    #length_row = length(poisd_sample_gene_of_interest)
    #bv_sum = sum(is.na(poisd_sample_gene_of_interest))
    #check_NA_matix = bv_sum - length_row
    #if (check_NA_matix>(-length_row)) {
    #  poisd_sample_gene_of_interest = 1- t(gene_list_of_interest) #parallelDist(t(gene_list_of_interest))
    #}

    #-----NOTE: It requires alot of memory ----!!!
    #poisd_sample_gene_corrected <- PoissonDistance(t(gene_exp_corrected))$dd
    #Check if this has NA
    #length_row = length(poisd_sample_gene_corrected)
    #bv_sum = sum(is.na(poisd_sample_gene_corrected))
    #check_NA_matix = bv_sum - length_row
    #if (check_NA_matix>(-length_row)) {
    #  poisd_sample_gene_corrected = 1 - cor(t(gene_exp_corrected))
    #  print("Sample distance: genes corrected: pearson distance ")
    #}else {
    #  print("Sample distance: genes corrected: poisson distance ")
    #}
    #---- Calculate the distance Euclidean -
    #dist_gene_of_interest = parallelDist(gene_list_of_interest)
    #dist_gene_corrected = parallelDist(gene_exp_corrected)

    ##---PAM ---
    #print("Perform: PAM Clustering...")
    #----PAM clustering (genes)-
    #print("First with genes...")
    #pamy_r_gene_of_interest <- pam(dist_gene_of_interest, 7)
    #kmrow_gene_of_interest <- pamy_r_gene_of_interest$clustering

    #pamy_r_gene_exp_corrected <- pam(dist_gene_corrected, 7)
    #kmrow_gene_exp_corrected <- pamy_r_gene_exp_corrected$clustering

    #---- PAM clustering (samples)-
    #print("Second with samples...")
    #pamy_c_gene_of_interest <- pam(as.dist(poisd_sample_gene_of_interest), 7)
    #kmcol_gene_of_interest <- pamy_c_gene_of_interest$clustering

    #pamy_c_gene_exp_corrected <- pam(as.matrix(poisd_sample_gene_corrected), 7)
    #kmcol_gene_exp_corrected <- pamy_c_gene_exp_corrected$clustering

    #---- Generate the dendrogram ----#
    #print("Generate a heat map based on poisson distance martix for samples and genese of selected genes...")
    #gene_list_q =rownames(gene_list_of_interest)
    #Gene/Sample (Selected)
    #hr_selected = hclust(as.dist(dist_gene_of_interest))
    #hc_selected = hclust(as.dist(poisd_sample_gene_of_interest))
    #hr_selected$labels = as.character(sampledata$SampleID)
    #hc_selected$labels = as.character(sampledata$SampleID)
    #Gene/Sample (All)
    #hr_all = hclust(as.dist(dist_gene_corrected))
    #hc_all = hclust(as.dist(poisd_sample_gene_corrected))
    #hr_all$labels = as.character(sampledata$SampleID)
    #hc_all$labels = as.character(sampledata$SampleID)
    #Genere unrooted tree
    #--- Selected genes --#
#    png(paste0(prefix, x_type,"_gene_selected_sample_hc_unrooted_tree.png"), width=20,height=20,units="in", res = 600)
#    plot(as.phylo(hc_selected), type = "unrooted",cex = 0.1)
#    dev.off()
#    png(paste0(prefix, x_type,"_gene_selected_gene_hc_unrooted_tree.png"), width=50,height=50,units="in", res = 600)
#    plot(as.phylo(hr_selected), type = "unrooted",cex = 0.1)
#    dev.off()
#    #--- all genes --#
#    png(paste0(prefix, x_type,"_gene_all_sample_hc_unrooted_tree.png"), width=20,height=20,units="in", res = 600)
#    plot(as.phylo(hc_all), type = "unrooted",cex = 0.1)
#    dev.off()
#    png(paste0(prefix, x_type,"_gene_all_gene_hc_unrooted_tree.png"), width=50,height=50,units="in", res = 600)
#    plot(as.phylo(hr_all), type = "unrooted",cex = 0.1)
#    dev.off()
####
    ##---Heatmap's presentation ---######
    Cluster_data = list()
    #selected genes
    Cluster_data$selected$name = "selected"
    Cluster_data$all$name = "all"
    ### ------------ Deactivated-------
    #Cluster_data$selected$hr = hr_selected ###
    #Cluster_data$selected$dist_row = dist_gene_of_interest
    #Cluster_data$selected$hc = hc_selected
    #Cluster_data$selected$dist_col = poisd_sample_gene_of_interest$dd
    #Cluster_data$selected$PAM_row = kmrow_gene_of_interest
    #Cluster_data$selected$PAM_col = kmcol_gene_of_interest
    #All genes

    #Cluster_data$all$hr = hr_all
    #Cluster_data$all$dist_row = dist_gene_corrected
    #Cluster_data$all$hc = hc_all
    #Cluster_data$all$dist_col = poisd_sample_gene_corrected$dd
    #Cluster_data$all$PAM_row = kmrow_gene_exp_corrected
    #Cluster_data$all$PAM_col = kmcol_gene_exp_corrected
    ###--------------Heatmap section -------------------######
    print(paste0("Number of all genes: ", nrow(gene_exp_corrected)))
    print(paste0("Number of selected genes: ", nrow(gene_list_of_interest)))


    #Data prepration
    df_class = data.frame(SampleType=as.factor(sampledata$sample_type), Batch = as.factor(sampledata$batch2), Batch_by_day= as.factor(sampledata$batch))
    rownames(df_class) = sampledata$SampleID

    ann_colors = list(
      SampleType = greenscale(length(levels(df_class$SampleType))), # c(T = "#7570B3", N = "#E7298A"),
      #SampleClass = redscale(max(x_itm$PAM_col)),
      Batch = blues9[1:length(levels(df_class$Batch))],
      Batch_by_day = rainbow(length(levels(df_class$Batch_by_day))),
      gene_class = rainbow(7))

    names(ann_colors$SampleType) = levels(df_class$SampleType)
    #names(ann_colors$SampleClass) = levels(df_class$SampleClass)
    names(ann_colors$Batch) = levels(df_class$Batch)
    names(ann_colors$Batch_by_day) = levels(df_class$Batch_by_day)

    for (x_itm in Cluster_data){
        ##-------- FOR ALL -----------####
        print(paste0("Proc: ", x_itm$name))
        #---Prepare the annotation data ---
        gene_for_heatmap = NULL
        if (x_itm$name=="selected") {
          gene_for_heatmap = gene_list_of_interest
        } else {
          gene_for_heatmap = gene_exp_corrected
        }
        print(dim(gene_for_heatmap))
        ## ---- Run COCA -----#
        print("Proc: Run COCA to cluster the samples")
        title=tempdir()
        results_COCA_R = ConsensusClusterPlus(gene_for_heatmap,maxK=20,reps=1000,pItem=0.8,pFeature=1, title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")
        icl = calcICL(results_COCA_R,title=title,plot="png")

        #Get the recommended cluster number
        recommend_cluster = get_cluster_number(icl=icl, prefix_id = prefix, normalization_id = x_type, data_type_id = x_itm$name)
        df_class$COCA = as.factor(results_COCA_R[[recommend_cluster]]$consensusClass)
        ann_colors$COCA = redgreen(recommend_cluster)
        cutree_col =results_COCA_R[[recommend_cluster]]$consensusTree
        cutree_col$labels = colnames(gene_for_heatmap)
        names(ann_colors$COCA) = levels(df_class$COCA)

        recommend_cluster_x = get_cluster_number_based_on_Consensus(icl=icl, prefix_id = prefix, normalization_id = x_type, data_type_id = x_itm$name)
        df_class$COCA_2 = as.factor(results_COCA_R[[recommend_cluster_x]]$consensusClass)
        ann_colors$COCA_2 = redgreen(recommend_cluster_x)
        names(ann_colors$COCA_2) = levels(df_class$COCA_2)

        #if (x_itm$name=="selected") {
        #  Cluster_data$selected$COCA =results_COCA_R
        #} else {
        #  Cluster_data$all$COCA =results_COCA_R
        #}

        #--- cluster gene using NbClust---#
        library("NbClust")
        set.seed(123)
        print("Perform Nbclust to determine the best partition of genes...(silhouette) (performing all 30 tests takes minutes till hours)")
        nb <- NbClust(gene_for_heatmap, distance = "euclidean", min.nc = 2,
                      max.nc = 20, method = "kmeans", index ="silhouette")
        #print(fviz_nbclust(nb) + theme_minimal())
        print(nb)
        ggsave(paste(prefix,x_type,x_itm$name,"_clustering.png",sep="_"))
        best_cluster = nb$Best.partition
        df_row = data.frame(gene_class = as.factor(best_cluster))
        rownames(df_row) = names(best_cluster)
        ann_colors$gene_class =  rainbow(length(levels(df_row$gene_class)))
        names(ann_colors$gene_class) = levels(df_row$gene_class)
        #detach("package:NbClust", unload=TRUE)
        
        #----cluster gene using tight cluster ---#
        if (x_itm$name=="selected") {
        print("Running: Tight cluster for gene clustering...")
        library(tightClust)

        gene_for_heatmap_c = gene_for_heatmap - rowMeans(gene_for_heatmap)
        gene_for_heatmap_c[gene_for_heatmap_c>4] <- 4
        gene_for_heatmap_c[gene_for_heatmap_c<(-4)] <- -4
        h_tightclust <- getTightCluster(x=gene_for_heatmap_c, standardize.gene = F)
        col_anno = data.frame(SampleType = df_class$SampleType)
        rownames(col_anno) = rownames(df_class)
        pheatmap_tclust <- generate_pheatmap_from_tight_cluster(h_tightclust, standardize.gene = F, order.sample = "SampleType", col_anno = col_anno)
        save_pheatmap_pdf(pheatmap_tclust, paste0(prefix,"_", x_type, "_", x_itm$name ,"_noise_filtered_tight_clust_gene_heatmap.pdf"), width =25, height=25)

        df_row$gene_tclust = as.factor(h_tightclust$cluster)
        ann_colors$gene_tclust = rainbow(length(levels(df_row$gene_tclust)))
        names(ann_colors$gene_tclust) = levels(df_row$gene_tclust)
        #detach("package:tightClust", unload=TRUE)
        }
        ##---- Deactivated memory issue -----#####
        #if (x_itm$name=="selected") {
        #  Cluster_data$selected$Voting_gene_Cluster =nb
        #  Cluster_data$selected$htclust_gene = h_tightclust
        #} else {
        #  Cluster_data$all$Voting_gene_Cluster =nb
        #  Cluster_data$selected$htclust = h_tightclust
        #}
        
        
        ##---cluster samples using mclust EM ----####
        #Consums alot of memory. Minimum 128 GB for >17000 genes.
        if (x_itm$name=="selected") {
        library(mclust)
        gene_g_top =t(gene_for_heatmap)
        BIC <- mclustBIC(gene_g_top)

        png(paste0(prefix,"_", x_type, "_", x_itm$name ,"_mclust_BIC_gene_plot.png"))
        plot(BIC)
        dev.off()

        print(summary(BIC))
        mod1 <- Mclust(gene_g_top, x = BIC)

        tbl_sample_distribution = table(mod1$classification)
        print(table(mod1$classification))
        print(length(mod1$classification))
        print(prop.table(tbl_sample_distribution)*100)

        print("Running: ICL for EM model...")
        ICL = mclustICL(gene_g_top)
        print(summary(ICL))
        png(paste0(prefix,"_", x_type, "_", x_itm$name ,"_mclust_gene_plot_icl.png"))
        plot(ICL)
        dev.off()
        #Add the cluster EM to the annotation
        df_class$sample_class_EM = as.factor(mod1$classification)
        ann_colors$sample_class_EM =  rainbow(length(levels(df_class$sample_class_EM)))
        names(ann_colors$sample_class_EM) = levels(df_class$sample_class_EM)
        #Remove the library
        #detach("package:mclust", unload=TRUE)
        }

        ##---- Generate Heatmap -----##
        print("Generate the heat map....")
        gene_for_heatmap_c  <- gene_for_heatmap - rowMeans(gene_for_heatmap)
        gene_for_heatmap_c[gene_for_heatmap_c>5] <- 4.9
        gene_for_heatmap_c[gene_for_heatmap_c<(-5)] <- -4.9
        #Get the pearson distance
        par_dis<-as.dist(1-cor(gene_for_heatmap_c));

        #Save the annotation
        col_anno2 = col_anno
        col_anno2$sample_id = rownames(col_anno)
        df_row2 = df_row
        df_row2$gene_id = rownames(df_row)
        write.table(col_anno2,paste0(prefix, "_", x_type, "_", x_itm$name, "_genes_sample_cluster_annotation.csv"), row.names=F, col.names=T)
        write.table(df_row2,paste0(prefix, "_", x_type, "_", x_itm$name, "_genes_gene_cluster_annotation.csv"), rownames=F, colnames=T)

        ptl_heatmap = pheatmap(gene_for_heatmap_c,
                               clustering_distance_cols = par_dis,
                               clustering_distance_rows = parallelDist(gene_for_heatmap_c),
                               annotation_row = df_row,
                               annotation_col = df_class,
                               fontsize = 6.5,
                               fontsize_row=1.4,
                               cluster_rows = 6,
                               fontsize_col = 6,
                               annotation_colors = ann_colors,
                               col = rev(brewer.pal(n = 11, name = "RdBu")))

        save_pheatmap_pdf(ptl_heatmap, paste0(prefix,"_", x_type, "_", x_itm$name ,"_gene_heatmap.pdf"), width =40, height=40)

        
        
        ## ----find 10% top variable genes and generate a heatmap based on that -------------------------------------------------------
        topVarGenes <- head(order(rowVars(gene_for_heatmap), decreasing = TRUE), round(nrow(gene_for_heatmap)*0.10))
        print(paste0("The number of top genes have been selected: ", round(nrow(gene_for_heatmap)*0.10)))
        file_name_x = paste0(prefix, "_", x_type, "_", x_itm$name, "_top_variable_genes.csv")
        tblx_to_savae = as.data.frame(gene_for_heatmap[topVarGenes,])
        tblx_to_savae$Gene = rownames(gene_for_heatmap[topVarGenes,])
        write.table(tblx_to_savae, file = file_name_x, col.names=T, sep=",", row.names = F)
        print(paste0("Top variable genes by ", x_itm$name, " genes using the normalization approach " , x_type))
        print(rownames(gene_for_heatmap[topVarGenes,]))
        mat  <- gene_for_heatmap[topVarGenes,]
        mat  <- mat - rowMeans(mat)
        mat[mat>4] <- 4
        mat[mat<(-4)] <- -4
        anno <- df_class
        
        #COCA
        #1. Samplewise
        print("----Samplewise----")
        title=tempdir()
        results_COCA_R = ConsensusClusterPlus(mat,maxK=20,reps=1000,pItem=0.8,pFeature=1, title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")
        icl = calcICL(results_COCA_R,title=title,plot="png")
        #Get the K number that has the lowest value.
        recommend_cluster = get_cluster_number(icl=icl, prefix_id = prefix, normalization_id = x_type, data_type_id = paste0(x_itm$name, "top_genes_samplewise"))
        anno$COCA = as.factor(results_COCA_R[[recommend_cluster]]$consensusClass)
        ann_colors$COCA = redgreen(recommend_cluster)
        names(ann_colors$COCA) = levels(anno$COCA)

        #2. Genewise
        print("----Genewise----")
        title=tempdir()
        results_COCA_R = ConsensusClusterPlus(t(mat),maxK=20,reps=1000,pItem=0.8,pFeature=1, title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")
        icl = calcICL(results_COCA_R,title=title,plot="png")
        #Get the K number that has the highest consensus.
        recommend_cluster = get_cluster_number(icl=icl, prefix_id = prefix, normalization_id = x_type, data_type_id = paste0(x_itm$name, "top_genes_genewise"))
        df_row = data.frame(gene_COCA = as.factor(results_COCA_R[[recommend_cluster]]$consensusClass))
        ann_colors$gene_COCA = blueyellow(recommend_cluster)
        names(ann_colors$gene_COCA) = levels(df_row$gene_COCA)

        recommend_cluster_x = get_cluster_number_based_on_Consensus(icl=icl, prefix_id = prefix, normalization_id = x_type, data_type_id = paste0(x_itm$name, "top_genes_genewise"))
        df_row$gene_COCA_2 = as.factor(results_COCA_R[[recommend_cluster_x]]$consensusClass)
        ann_colors$gene_COCA_2 = blueyellow(recommend_cluster_x)
        names(ann_colors$gene_COCA_2) = levels(df_row$gene_COCA_2)

        #3. Gene 2. approach nbclust
        set.seed(123)
        library(NbClust)
        print("Perform Nbclust to determine the best partition of genes (silhouette)... (Takes minutes till hours)")
        nb <- NbClust(mat, distance = "euclidean", min.nc = 2,
                      max.nc = 20, method = "kmeans", index ="silhouette")
        print(nb)
        #print(fviz_nbclust(nb) + theme_minimal())
        ggsave(paste(prefix,x_type,x_itm$name, "_top_gene_clustering.png",sep="_"))
        #library("fpc")
        best_cluster = nb$Best.partition
        df_row$gene_class = as.factor(best_cluster)
        rownames(df_row) = names(best_cluster)
        ann_colors$gene_class =  rainbow(length(levels(df_row$gene_class)))
        names(ann_colors$gene_class) = levels(df_row$gene_class)
        #detach("package:NbClust", unload=TRUE)
        
        
        #3. Gene 3. approach tight cluster
        #if (x_itm$name=="selected") {
        library(tightClust)
        print("Running: Tight cluster for gene clustering...")
        h_tightclust <- getTightCluster(x=mat)
        col_anno = data.frame(SampleType = anno$SampleType)
        rownames(col_anno) = rownames(df_class)
        pheatmap_tclust <- generate_pheatmap_from_tight_cluster(h_tightclust, standardize.gene = F, order.sample = "SampleType", col_anno = col_anno)
        save_pheatmap_pdf(pheatmap_tclust, paste0(prefix,"_", x_type, "_", x_itm$name ,"_noise_filtered_tight_clust_top_gene_heatmap.pdf"), width =10, height=10)
        length(h_tightclust$cluster)
        df_row$gene_tclust = as.factor(h_tightclust$cluster)
        ann_colors$gene_tclust = rainbow(length(levels(df_row$gene_tclust)))
        names(ann_colors$gene_tclust) = levels(df_row$gene_tclust)
        #detach("package:tightClust", unload=TRUE)
        #}
        #Deactivate because of memory issue.
        #Save the results.
        #if (x_itm$name=="selected") {
        #  Cluster_data$selected$top_gene_voting_gene_Cluster =nb
        #  Cluster_data$selected$top_gene_voting_gene_tightCluster =h_tightclust
        #} else {
        #  Cluster_data$all$top_gene_voting_gene_Cluster =nb
        #  Cluster_data$all$top_gene_voting_gene_Cluster = pheatmap_tclust
        #}

        #if (x_itm$name=="selected") {
        #4. Gene 4. model-based clustering (Gaussian finite mixture model fitted by EM algorithm)
        print("Running: Model-based clustering using mclust...")
        library(mclust)
        gene_g_top = t(mat)
        BIC <- mclustBIC(gene_g_top)

        png(paste0(prefix,"_", x_type, "_", x_itm$name ,"_mclust_BIC_top_gene_plot.png"))
        plot(BIC)
        dev.off()

        print(summary(BIC))
        mod1 <- Mclust(gene_g_top, x = BIC)
        tbl_sample_distribution = table(mod1$classification)
        print(table(mod1$classification))
        print(prop.table(tbl_sample_distribution)*100)

        print("Running: ICL for EM model...")
        ICL = mclustICL(gene_g_top)
        print(summary(ICL))
        png(paste0(prefix,"_", x_type, "_", x_itm$name ,"_mclust_top_gene_plot_icl.png"))
        plot(ICL)
        dev.off()
        #Add the cluster EM to the annotation
        anno$sample_class_EM = as.factor(mod1$classification)
        ann_colors$sample_class_EM =  rainbow(length(levels(anno$sample_class_EM)))
        names(ann_colors$sample_class_EM) = levels(anno$sample_class_EM)
        #detach("package:mclust", unload=TRUE)
        #}

        #5. Save the annotation samples and gene
        col_anno2 = anno
        col_anno2$sample_id = rownames(anno)
        df_row2 = df_row
        df_row2$gene_id = rownames(df_row)
        write.table(anno,paste0(prefix, "_", x_type, "_", x_itm$name, "_top_genes_sample_cluster_annotation.csv"), sep=",", row.names=F, col.names=T)
        write.table(df_row,paste0(prefix, "_", x_type, "_", x_itm$name, "_top_genes_gene_cluster_annotation.csv"), sep=",", rownames=F, colnames=T)
        #6. Generate the heatmap
        print("Generating heatmap for top variable genes....")
        ptl_heatmap = pheatmap(mat, annotation_col = anno, annotation_colors= ann_colors,
                               annotation_row=df_row, clustering_distance_rows="euclidean",
                               clustering_distance_cols="correlation",
                               fontsize = 6.5,
                               fontsize_row=3,
                               fontsize_col = 2,
                               col = rev(brewer.pal(n = 11, name = "RdBu")))
        save_pheatmap_pdf(ptl_heatmap, paste0(prefix, "_", x_type, "_", x_itm$name,"_top_variable_genes_heatmap.pdf"), width=10,height=10)
    }
  }
  print("Done...")
  save.image(file =  paste0(prefix, "_results.RData"))
  print("Save the memory image to the hard disk...")

  ## --- What is the next step? ----
  if (PlayWithResults==T) {
    print("Run: Further analyses, heatmap....")
    Y$Report <- RunStandardResultsReport(Y, contrast_list=contrast_list,group_factor=group_factor, prefix = prefix)
    return(Y)
  } else {
  return(Y) }
}
### ---- RunStandardResultsReport -----
RunStandardResultsReport <- function(result, contrast_list=c("sample_type", "T", "N"), group_factor=c("sample_type","batch", "SampleID"), int_group= c("sample_type","batch"), prefix="") {
  print(paste0("Run: ", contrast_list))
  for (itm_c in contrast_list) {
    for (itm_b in contrast_list) {
      if (itm_c != itm_b) {
        dds <- result$dds
        con_group_list = c("sample_type", itm_b, itm_c)
        con_groupt_txt = paste0(itm_b, "_v_", itm_c)
        ##--- Save the results of the selected subgroup ----
        print("Run: prepare DES results...")
        res <- results(result$dds, contrast=con_group_list)
        print("Save: prepare DES results")
        write.csv(as.data.frame(res),
        file=paste0(prefix,"_", con_groupt_txt,"_resDES_results.csv"))
        ##--- Run the independent hypothesis weighting and save the result... ----
        print("Run: Independent Hypothesis Weighting")
        resIHW <- results(result$dds, filterFun=ihw)
        print(summary(resIHW))
        print(sum(resIHW$padj < 0.1, na.rm=TRUE))
        print("Save: Independent Hypothesis Weighting")
        write.csv(as.data.frame(resIHW), file=paste0(prefix,"_", con_groupt_txt,"resIHW_condition_results.csv"))
        print("Done: Independent Hypothesis Weighting")
        ## ---- Effects of transformations on the variance ----
        print(paste0("Run: ", "Effects of transformations on the variance"))
        notAllZero <- (rowSums(counts(dds))>0)
        png(paste0(prefix, "_", con_groupt_txt, "_plot_log2.png"))
        meanSdPlot(log2(counts(dds,normalized=TRUE)[notAllZero,] + 1))
        dev.off

        ## ----- Heatmap of the count matrix --------
        select <- order(rowMeans(counts(dds,normalized=TRUE)),
                    decreasing=TRUE)
        nt <- normTransform(dds) # defaults to log2(x+1)
        log2.norm.counts <- assay(nt)[select,]
        df <- as.data.frame(colData(dds)[,group_factor])
        ptl_heatmap = pheatmap(log2.norm.counts, cluster_rows=FALSE, show_rownames=FALSE,cluster_cols=FALSE, annotation_col=df)
        save_pheatmap_pdf(ptl_heatmap, paste0(prefix,"_", con_groupt_txt, "_gene_log2_normalized_heatmap.pdf"))
        ## ------------ Run mcols. ------------------------------------------
        print("Run mcols.")
        print(mcols(res, use.names = TRUE))
        print("Print a summary...")
        print(summary(res))

        print("alpha =0.05 and padj<0.05")
        res.05 <- results(dds, alpha = 0.05)
        print(table(res.05$padj < 0.05))

        print("resLFC1 <- results(dds, lfcThreshold=1) and padj<0.1")
        resLFC1 <- results(dds, lfcThreshold=1)
        print(table(resLFC1$padj < 0.1))

        print("sum(res$pvalue < 0.05, na.rm=TRUE)")
        print(sum(res$pvalue < 0.05, na.rm=TRUE))
        print("sum(!is.na(res$pvalue))")
        print(sum(!is.na(res$pvalue)))
        print("sum(res$padj < 0.1, na.rm=TRUE)")
        print(sum(res$padj < 0.1, na.rm=TRUE))

        resSig <- subset(res, padj < 0.1)
        print("head(resSig[ order(resSig$log2FoldChange), ]) & padj<0.01")
        print(head(resSig[ order(resSig$log2FoldChange), ]))

        print("head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ])")
        print(head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ]))

        ## ----plotcounts----------------------------------------------------------
        topGene <- rownames(res)[which.min(res$padj)]
        png(paste0(prefix, "_", con_groupt_txt, "_plot_counts.png"))
        plotCounts(dds, gene = topGene, intgroup=int_group)
        dev.off()
        ## ----ggplotcountsjitter, fig.width = 4, fig.height = 3-------------------
        geneCounts <- plotCounts(dds, gene = topGene, intgroup = int_group,
                           returnData = TRUE)
        print(colnames(geneCounts))
        ggplot(geneCounts, aes(x = sample_type, y = count, color = batch)) +
        scale_y_log10() +  geom_beeswarm(cex = 3)
        ggsave(paste0(prefix,"_", con_groupt_txt,"_geneplotCounts_sample_type_batch.png"))
        ## ----ggplotcountsgroup, fig.width = 4, fig.height = 3--------------------
        ggplot(geneCounts, aes(x = sample_type, y = count, color = batch, group = batch)) +
        scale_y_log10() + geom_point(size = 3) + geom_line()
        ggsave(paste0(prefix,"_", con_groupt_txt,"_geneplotCounts_sample_type_batch.png"))
        ## ----plotmaNoShr---------------------------------------------------------
        res.noshr <- results(dds)
        png(paste0(prefix, "_", con_groupt_txt,"_plotma_no_shr.png"))
        plotMA(res.noshr, ylim = c(-5, 5))
        dev.off()
        ## ----plotma--------------------------------------------------------------
        res <- lfcShrink(dds, contrast=con_group_list, res=res)
        print(res)
        png(paste0(prefix,"_", con_groupt_txt,"_plotma.png"))
        plotMA(res, ylim = c(-5, 5))
        dev.off()
        ## ----plotmalabel---------------------------------------------------------
        png(paste0(prefix, "_", con_groupt_txt,"_plotma_label.png"))
        plotMA(res, ylim = c(-5,5))
        topGene <- rownames(res)[which.min(res$padj)]
        with(res[topGene, ], {
          points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
          text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
        })
        dev.off()
        ## ----histpvalue2---------------------------------------------------------
        png(paste0(prefix, "_", con_groupt_txt,"_histogram.png"))
        hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
             col = "grey50", border = "white")
        dev.off()
        ## ----sensitivityovermean, fig.width=6------------------------------------
        png(paste0(prefix,"_", con_groupt_txt, "_sensitivityovermean.png"))
        qs <- c(0, quantile(resLFC1$baseMean[resLFC1$baseMean > 0], 0:6/6))
        bins <- cut(resLFC1$baseMean, qs)
        levels(bins) <- paste0("~", round(signif((qs[-1] + qs[-length(qs)])/2, 2)))
        fractionSig <- tapply(resLFC1$pvalue, bins, function(p)
        mean(p < .05, na.rm = TRUE))
        barplot(fractionSig, xlab = "mean normalized count",
          ylab = "fraction of small p values")
        dev.off()

        library("ReportingTools")
        htmlRep <- HTMLReport(shortName=prefix, title=prefix,
                        reportDirectory= paste0("./", prefix,"_", con_groupt_txt,"_report"))
        res_d_d = as.data.frame(res)
        res_d_d$gene = rownames(res)
        publish(as.data.frame(res_d_d), htmlRep)
        url <- finish(htmlRep)
      }
    }
  }
  print("Save the memory image to the hard disk...")
  save.image(file =  paste0(prefix, "_results_presentation.RData"))
  return(result)
}
