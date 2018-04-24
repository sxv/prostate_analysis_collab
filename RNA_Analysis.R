## try http:// if https:// URLs are not supported
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
    biocLite(new.pkg)
  sapply(pkg, require, character.only = TRUE)
}

if("slam" %in% rownames(installed.packages()) == FALSE) {
  install.packages('devtools')
  library(devtools)
  slam_url <- "https://cran.r-project.org/src/contrib/Archive/slam/slam_0.1-37.tar.gz"
  install_url(slam_url)}

packages<-c("som","dplyr", "ggbeeswarm", "readr","ggplot2", "pheatmap", "magrittr", "MASS", "psych", "RColorBrewer", "ape")
biopackages <- c("preprocessCore", "BiocParallel", "PoiClaClu", "ConsensusClusterPlus", "GenomicFeatures", "DESeq2", "vsn", "sva", "ReportingTools", "IHW","Glimma")
check.packages(packages)
check.biocLites(biopackages)

library(preprocessCore)#
#library(som)#
library(dplyr)#
library(ConsensusClusterPlus)#
library(ggplot2)#
#library(MASS) #
library(PoiClaClu)#
library("pheatmap")#
library("RColorBrewer")#
library("GenomicFeatures")#
library("DESeq2")#
library(psych)#
#library(sva)#
library("magrittr")#
#library(vsn)#
library("BiocParallel")#
library("IHW")
library(readr)
library(ape)
#####
# Set hyperparameter
rm(list=ls())
initial_filter = "rowSum" # rowmean rowSum none (Standard)
cutoff_level = 1
select_only_exon_genes = T
#####
#functions
#####
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
#####
#Set a work directory and files
#####
#dir_path = "/Users/okyazeminaga/OneDrive/Forschung/RNA_Seq/"
#setwd(dir_path)
rna_seq_data_filename = "prostate225.reads.sorted_patient.simple.csv"
patient_sample_annotation_filename = "prostate225.design_file.batched_patient.simple.csv"
#######
# Data preparation
######
#Load data
data_seq_raw_a = read.csv(rna_seq_data_filename)
sample_data <- read.table(patient_sample_annotation_filename, sep=",", header=T, stringsAsFactors = F)
###
#Define samples were separately sequenced at different times.
list_batch1 = c("32","33","34","35","36","37")
list_batch2 = c("N","P","S","T","U","V","W", "X")
###
#Reshape the matrix by moving the gene name to the row name
data_seq_raw = data_seq_raw_a[,-1]
row.names(data_seq_raw) = data_seq_raw_a[,1]
data_seq_raw_matrix = data.matrix(data_seq_raw)
print("Show the data dimension (genexsample):")
print(paste0(ncol(data_seq_raw_matrix), "x", nrow(data_seq_raw_matrix)))
print("Show the dimension of the sample annotation")
print(paste0(nrow(sample_data), "x", ncol(sample_data)))
#####
#Change the sample annotation data
#####
sample_data$id[sample_data$id=="9696.35.T43"] = "9696.35.T34"
sample_data$id[sample_data$id=="W88.NA"] = "W88.NN"
sample_data$id[sample_data$id=="W89.NA" ] = "W89.NN"
sample_data$id[sample_data$id=="W90.NA"] = "W90.NN"
sample_data$id[sample_data$id=="P109.NA"] = "P109.NN"
#####
# Get the sample id from gene data
#####
#Case labeling by pathologic finding
colnames(data_seq_raw_matrix) = gsub("NA", "NN", colnames(data_seq_raw_matrix))
colnames(data_seq_raw_matrix) = gsub(".T43", ".T34", colnames(data_seq_raw_matrix))
dt = strsplit(colnames(data_seq_raw_matrix), '\\.')
name_ls = list()
case_ls = list()
for (i in 1:length(dt)) {
  v = unlist(dt[i])
  len_c = length(v)
  name_ls[i] = v[len_c]
  p_n = ""
  for (pl in 1:(len_c-1)) {
    if (!(grepl("X_", v[pl]))) {
      if (tolower(v[pl])!= "normal"){
        if (pl >1 & nchar(p_n)>0) {
          p_n = paste0(p_n, ".")
        }

        p_n = paste0(p_n, v[pl])
      }
    }
  }
  if (grepl("V", p_n)) {
    p_n = "V"
  }
  if (grepl("N", p_n)) {
    p_n = "N"
  }
  if (grepl("S", p_n)) {
    p_n = "S"
  }
  if (grepl("T", p_n)) {
    p_n = "T"
  }
  if (grepl("U", p_n)) {
    p_n = "U"
  }
  if (grepl("W", p_n)) {
    p_n = "W"
  }
  if (grepl("P", p_n)) {
    p_n = "P"
  }
  if (grepl("X",p_n)) {
    p_n = "X"
  }
  case_ls[i] = p_n
}

sample_types = unlist(name_ls)
print(paste0("Sample number: ", length(sample_types)))
case_assigment = unlist(case_ls)
case_id = unique(case_assigment)
print(paste0("The number of cases: ", length(case_id)))
#Identify unique findings
q = unique(sample_types)
tbl_sample_distribution = table(sample_types)
tbl_sample_distribution
prop.table(tbl_sample_distribution)*100
#Check if these groups have well distributed number of samples -Ignore-
ks.test(as.numeric(as.factor(sample_types)),"pnorm")
ks.test(as.numeric(as.factor(sample_types)),"pgamma", 3, 2)
ks.test(as.numeric(as.factor(sample_types)),"dpois", 3, 2)
#Modify the sample id to make consistent cross the data
#Add possible batch effects because of different sequencing times
sample_data$batch2 = "0"
sample_data$batch2[sample_data$case %in% list_batch2] = "1"
sample_data$sample_types = sample_types
#####
# Filter and Ploting the row distribution
#####
#Run the initial filter function
print("Run the initial filter...")
data_seq_raw_matrix_reduced = NULL
if (tolower(initial_filter)=="rowmean") {
  data_seq_raw_matrix_reduced = data_seq_raw_matrix[rowMeans(data_seq_raw_matrix) >=cutoff_level, ]
  print(paste0("gene number after the initial filter: ", nrow(data_seq_raw_matrix_reduced)))
  print(paste0("reduced by: ", nrow(data_seq_raw_matrix)-nrow(data_seq_raw_matrix_reduced)))
}

if (is.null(data_seq_raw_matrix_reduced)) {
  data_seq_raw_matrix_reduced = data_seq_raw_matrix
}
#illustrate the data distribution--> likly poisson distribution and gamma distribution, chiq-square distribution
#illustrate the data distribution--> likly poisson distribution and gamma distribution, chiq-square distribution
#data_seq_raw_matrix = data_seq_raw_matrix[rowMeans(data_seq_raw_matrix) >=1, ]
colramp = colorRampPalette(c(3,"white",2))(ncol(data_seq_raw_matrix_reduced))
#Distribution analyses
plot(density(data_seq_raw_matrix_reduced[,1]),col=colramp[1],lwd=3,ylim=c(0,4), xlim=c(0,1000))
for(i in 2:ncol(data_seq_raw_matrix_reduced)){lines(density(data_seq_raw_matrix_reduced[,i]),lwd=3,col=colramp[i])}
lines(density(rowMeans(data_seq_raw_matrix_reduced)), lwd=3, col="blue")
lines(density(apply(data_seq_raw_matrix_reduced, 1, median)), lwd=3, col="Yellow")
lines(density(apply(data_seq_raw_matrix_reduced, 1, geometric.mean)), lwd=3, col="black")
#Cumulative Distribution function
colramp = colorRampPalette(c(3,"white",2))(ncol(data_seq_raw_matrix_reduced))
P = ecdf(data_seq_raw_matrix_reduced)
plot(P, main="Cumulative Distribution function", xlab="raw read count", ylab="percentile")
dim(data_seq_raw_matrix_reduced)
######
# Run DESEQ2
######
#Prepare the factors (1. All 2. Binary: Normal vs. Tumor, 3. Trinty: Normal, Precursor, Tumor)
sample_data$sample_type_factors = factor(sample_data$sample_types, levels=c("NN", "NervN", "ST", "BPH", "LGPIN", "HGPIN", "T3", "T34", "T4", "T45", "T5", "TD", "PNI", "SVI", "LNM"))
sample_data$sample_type_boolean_factors = "T"
sample_data$sample_type_boolean_factors[sample_data$sample_types %in% c("NN", "NervN", "ST", "BPH", "LGPIN", "HGPIN")] = "N"
sample_data$sample_type_boolean_factors = factor(sample_data$sample_type_boolean_factors, levels=c("N", "T"))
sample_data$sample_type_triple_factors = "T"
sample_data$sample_type_triple_factors[sample_data$sample_types %in% c("LGPIN", "HGPIN")] = "P"
sample_data$sample_type_triple_factors[sample_data$sample_types %in% c("NN", "NervN", "ST", "BPH")] = "N"
sample_data$sample_type_triple_factors = factor(sample_data$sample_type_triple_factors, levels=c("N", "P" ,"T"))
#####
# Overall analyses
#####
# Select only exon genes
#####
if (select_only_exon_genes) {
  gene_data <- read_csv("hg19v2_clean.csv")
  gene_selected <-gene_data$gene_name[gene_data$biotype=="protein_coding"] #read.csv2("gene_list.csv", header = T)
  gene_data_raw_for_DESEQ2 = data_seq_raw_matrix_reduced[rownames(data_seq_raw_matrix_reduced) %in% gene_selected,]
} else {
  gene_data_raw_for_DESEQ2 = data_seq_raw_matrix_reduced
}
gene_data_raw_for_DESEQ2 = na.omit(gene_data_raw_for_DESEQ2)

dim(gene_data_raw_for_DESEQ2)
#####
#Annotation data for samples
#1. Binary factor
#Mandatory: A dataframe should have "sample_type" (as.factor) and "batch" before using it for RunSEQ
#USE data_seq_raw_matrix_reduced to retrieve the expression data
factors_to_run = c("sample_type_boolean_factors", "sample_type_triple_factors", "sample_type_factors")
factor_X = "sample_type_boolean_factors"
counter_id=1
#General move
source("RunSEQ.R")
factor_X = "sample_type_boolean_factors"
database_resutls_list = list()
for (factor_X in factors_to_run) {
  sample_Data_for_seq = data.frame(SampleID=colnames(gene_data_raw_for_DESEQ2), sample_type=sample_data[factor_X], batch = sample_data$batch2, batch2 = sample_data$batch, caseid = sample_data$case)
  colnames(sample_Data_for_seq)[colnames(sample_Data_for_seq)==factor_X] = "sample_type"
  database_resutls_list[counter_id] = RunSEQ(x=gene_data_raw_for_DESEQ2, initial_filter = initial_filter, set_cutoff = cutoff_level, Col_Data=sample_Data_for_seq, prefix=factor_X, Run_Id=counter_id, description_txt=factor_X, PlayWithResults=T, contrast_list=c(levels(sample_Data_for_seq$sample_type)))
  #cat(factor_X, out, file="RNA_Seq_Exp_log.txt", sep="n", append=TRUE)
  counter_id = counter_id +1
  save.image(file =  paste0(factor_X, "_final_stage_results_presentation.RData"))
}
save.image(file = "database_resutls_list_final_stage_results_presentation.RData")
#Case/sample specific
to_investigate_sample = colnames(model.matrix(~ sample_type, sample_Data_for_seq))
to_investigate_sample = to_investigate_sample[-1]

#Comparsion analyses
list_gene_findings_upregulated_sample = list()
list_gene_findings_downregulated_sample = list()
list_gene_findings_sig_gene_expression = list()
sample_Data_for_seq = data.frame(SampleID=colnames(gene_data_raw_for_DESEQ2), sample_type=sample_data$sample_type_factors, batch = sample_data$batch2, batch2 = sample_data$batch, caseid = sample_data$case)

dds <- DESeqDataSetFromMatrix(countData = gene_data_raw_for_DESEQ2, colData = sample_Data_for_seq, design= ~ sample_type + sample_type:caseid)
result_x <- RunSEQ(dds=dds, initial_filter = initial_filter, set_cutoff = cutoff_level, Col_Data=sample_Data_for_seq, prefix=factor_X, Run_Id=counter_id, description_txt=factor_X, PlayWithResults=Y, contrast_list=NA)

#Independent filtering of results
png("dds_disp_estimate_x.png")
plotDispEsts(dds)
dev.off
png("dds_number of rejections_x.png")
plot(metadata(result_x)$filterNumRej,
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(resuls_x)$lo.fit, col="red")
abline(v=metadata(result_x)$filterTheta)
dev.off
for (itm_1 in to_investigate_sample) {
  for (itm_2 in to_investigate_sample) {
    if (itm_1 != itm_2) {
    value_crit = paste0(itm_1,"_", itm_2)
    res = results(result_x$dds, contrast=list(itm_1,itm_2))
    sum(res$padj < 0.1, na.rm=TRUE)
    resSig = resSig[ order(-resSig$log2FoldChange ),]
    #Upregulated genes
    print("genes signficantly upregulated...")
    resSig_upregulated <- subset(resSig, (log2FoldChange >= 1))
    resSig_upregulated = resSig_upregulated[ order(-resSig_upregulated$log2FoldChange ),]
    Upregulated_gene_list = row.names(resSig_upregulated)
    Upregulated_gene_list
    #Downregulated genes
    print("genes signficantly downregulated...")
    resSig_downregulated <- subset(resSig, (log2FoldChange <= -1))
    resSig_downregulated = resSig_downregulated[ order(resSig_downregulated$log2FoldChange ),]
    Downregulated_gene_list = row.names(resSig_downregulated)
    Downregulated_gene_list
    #Store the results of analyses
    list_gene_findings_upregulated[[comparsion_string]] = Upregulated_gene_list
    list_gene_findings_downregulated[[comparsion_string]] = Downregulated_gene_list
    list_gene_findings_sig_gene_expression[[comparsion_string]] = resSig

    #MA plot
    summary(res)
    png(paste0(value_crit, prefix_x,"_MA_plot.png"))
    ptl_MA =plotMA(res, ylim=c(min(res$log2FoldChange[is.na(res$log2FoldChange)==F]),max(res$log2FoldChange[is.na(res$log2FoldChange)==F])))
    ptl_MA
    dev.off()
    }
  }
}
#Up regulated
col_list = list()
lenghth_list = lapply(list_gene_findings_upregulated, function(x) length(x))
length_list = unlist(lenghth_list)
value_max = max(length_list) +5
for (i in 1:length(list_gene_findings_upregulated)) {
  col_list = c(col_list, names(list_gene_findings_upregulated[i][1]))

}
col_list = unlist(col_list)

counter = 1
string_x = list()
df = NA
for (itm in list_gene_findings_upregulated) {
  max_local = length(itm) +1
  itmx = itm
  for (i in max_local:value_max) {
    itmx = c(itmx, "")
  }
  if (is.na(df)) {
    df = data.frame(itmx)
    colnames(df) = col_list[counter]
  }
  else {
    df[[col_list[counter]]] <- itmx
  }
  counter = counter + 1
}
#View(df)
write.table(df, 'comparison_analyes_upregulated_genes.csv', col.names = TRUE, append= F, sep=',', row.names = F)

#Down regulated
col_list = list()
lenghth_list = lapply(list_gene_findings_downregulated, function(x) length(x))
length_list = unlist(lenghth_list)
value_max = max(length_list) +5
for (i in 1:length(list_gene_findings_downregulated)) {
  col_list = c(col_list, names(list_gene_findings_downregulated[i][1]))

}
col_list = unlist(col_list)

counter = 1
string_x = list()
df = NA
for (itm in list_gene_findings_downregulated) {
  max_local = length(itm) +1
  itmx = itm
  for (i in max_local:value_max) {
    itmx = c(itmx, "")
  }
  if (is.na(df)) {
    df = data.frame(itmx)
    colnames(df) = col_list[counter]
  }
  else {
    df[[col_list[counter]]] <- itmx
  }
  counter = counter + 1
}
#View(df)
write.table(df, 'comparison_analyes_downregulated_genes.csv', col.names = TRUE, append= F, sep=',', row.names = F)
#All
col_list = list()
for (i in 1:length(list_gene_findings_sig_gene_expression)) {
  col_list = c(col_list, names(list_gene_findings_sig_gene_expression[i][1]))
}
col_list = unlist(col_list)

counter = 1
string_x = list()
for (itm in list_gene_findings_sig_gene_expression) {
  df = data.frame(itm)
  #colnames(df) = c("Gene", colnames(df))
  #counter = counter + 1
  df = tibble::rownames_to_column(df, var = "Gene")
  write.table(df, paste0(col_list[counter], '_comparison_analyes_all_sig_genes.csv'), col.names = TRUE, row.names = F, append= F, sep=',' )
  counter = counter + 1
}
save.image(file =  paste0("After_comparison_results_presentation.RData"))

###
title=tempdir()
nt <- normTransform(dds)
log2.norm.counts <- assay(nt)
modcombat = model.matrix(~as.factor(sample_type), data=sample_Data_for_seq)
log2.norm.counts.corrected = ComBat(dat=log2.norm.counts, batch=batch, mod=modcombat, par.prior=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
vsd.norm.counts <- assay(vsd)
vsd.norm.counts.corrected = ComBat(dat=vsd.norm.counts, batch=batch, mod=modcombat, par.prior=TRUE)
ptl_heatmap = pheatmap(log2.norm.counts.corrected , cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
save_pheatmap_pdf(ptl_heatmap, paste0(prefix,"_gene_log2_norm_corrected_heatmap.pdf"))
ptl_heatmap = pheatmap(vsd.norm.counts.corrected, cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
save_pheatmap_pdf(ptl_heatmap, paste0(prefix,"_gene_svd_corrected_heatmap.pdf"))
#COCA
results_COCA_log2_norm = ConsensusClusterPlus(log2.norm.counts.corrected,maxK=20,reps=1000,pItem=0.8,pFeature=1, title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")
results_COCA_log2_vsd = ConsensusClusterPlus(vsd.norm.counts.corrected,maxK=20,reps=1000,pItem=0.8,pFeature=1, title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")
#Cluster
dd=as.dist(1- cor(log2.norm.counts.corrected))
hc = hclust(dd)
clus5 = cutree(hc, 24)

plot(as.phylo(hc), type = "unrooted",cex = 0.1,  tip.color=rainbow(24)[clus5])
# fan
plot(as.phylo(hc), type = "fan", cex = 0.3,  tip.color=rainbow(24)[clus5])
# cladogram
plot(as.phylo(hc), type = "cladogram", cex = 0.3,  tip.color=rainbow(24)[clus5])
# radial
plot(as.phylo(hc), type = "radial", cex= 0.2, tip.color=rainbow(24)[clus5])
#VSD
dd_VST=as.dist(1- cor(vsd.norm.counts.corrected))
hc_VST = hclust(dd_VST)
clus5_VST = cutree(hc_VST, 24)

plot(as.phylo(hc_VST), type = "unrooted",cex = 0.1,  tip.color=rainbow(24)[clus5_VST])
# fan
plot(as.phylo(hc_VST), type = "fan", cex = 0.3,  tip.color=rainbow(24)[clus5_VST])
# cladogram
plot(as.phylo(hc_VST), type = "cladogram", cex = 0.3,  tip.color=rainbow(24)[clus5_VST])
# radial
plot(as.phylo(hc_VST), type = "radial", cex= 0.2, tip.color=rainbow(24)[clus5_VST])

## ------------------------------------------------------------------------
devtools::session_info()
save.image(file =  paste0(prefix, "final_stage_results_presentation.RData"))
##
