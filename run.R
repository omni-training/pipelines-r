#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(jsonlite)
  library(anndataR)
})

parser <- ArgumentParser(description = "Benchmarking entrypoint")

# dummy args (old API)
parser$add_argument("--output_dir", dest = "output_dir")
parser$add_argument("--name", dest = "name")

# inputs
parser$add_argument("--data_h5ad", dest = "data_h5ad", type = "character",
                    help = "input data h5ad path", required = TRUE)
parser$add_argument("--data_clusters_truth_num", dest = "data_clusters_truth_num", 
                    type = "character", help = "input true number of clusters path", required = TRUE)

# params
parser$add_argument("--method_name", dest = "method_name", type = "character",
                    help = "name of the method", choices = c("osca", "scrapper", "seurat", "bixverse"), required = TRUE)
parser$add_argument("--d_cluster", dest = "d_cluster", type = "integer",
                    help = "delta to true number of clusters", required = TRUE)
parser$add_argument("--n_comp", dest = "n_comp", type = "integer",
                    help = "number of PCA components to use for KNN graph construction",
                    default = 50, required = FALSE)
parser$add_argument("--n_neig", dest = "n_neig", type = "integer",
                    help = "number of neighbors to use for KNN graph construction", default = 15, required = FALSE)
parser$add_argument("--n_hvg", dest = "n_hvg", type = "integer",
                    help = "number of highly variable genes to use", default = 1000, required = FALSE)
parser$add_argument("--filter", dest = "filter", type = "character",
                    help = "filtering (manual uses suggested cutoffs, auto uses package QC)",
                    choices = c("manual", "auto"), default = "manual", required = FALSE)

# outputs
parser$add_argument("--output", action = "append", default = list(),
                    help = "Repeatable: key=path")
args <- parser$parse_args()
output <- setNames(sub("^[^=]*=", "", args$output),
                   sub("=.*$", "", args$output))

# logging
cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))

cat(sprintf("LOG: command line args\n----------------------------------\n"))
for (i in 1:length(args)) {
  cat(sprintf("  %s: %s\n", names(args)[i], args[[i]]))
}
for (i in 1:length(output)) {
  cat(sprintf("[output] %s: %s\n", names(output)[i], output[[i]]))
}
cat(sprintf("----------------------------------\n"))


# pickup helper code
cargs <- commandArgs(trailingOnly = FALSE)
m <- grep("--file=", cargs)
run_dir <- dirname(gsub("--file=", "", cargs[[m]]))
search_res_path <- file.path(run_dir, "search_res.r")
source(search_res_path)


# read data in
sce <- read_h5ad(args$data_h5ad, as = "SingleCellExperiment")
names(assays(sce))[3] <- "counts"
clusters_truth_num <- as.integer(readLines(args$data_clusters_truth_num, n = 1))
n_cluster <- clusters_truth_num + args$d_cluster

# define slots for outputs
time <- list(gpu_load = NA_real_, find_mit_gene = NA_real_, 
             filter = NA_real_, normalization = NA_real_, hvg = NA_real_, 
             scaling = NA_real_, pca = NA_real_, t_sne = NA_real_, 
             umap = NA_real_, louvain = NA_real_, leiden = NA_real_)

clustering_info <- list(resolutions = list(louvain = NA_real_, leiden = NA_real_), 
                                           num_runs = list(louvain = NA_integer_, 
                                                           leiden = NA_integer_))

# source and run appropriate method
if (args$method_name == "seurat") {
  # seurat ----------------------------------------------------------------
  source(file.path(run_dir, "seurat.R"))
  output_data <- run_seurat(sce, n_cluster, args$n_comp, args$n_neig, 
                            args$n_hvg, args$filter, time, clustering_info)
} else if (args$method_name == "osca") {
  # OSCA ------------------------------------------------------------------
  source(file.path(run_dir, "OSCA.R"))
  output_data <- run_osca(sce, n_cluster, args$n_comp, args$n_neig, 
                          args$n_hvg, args$filter, time, clustering_info)
} else if (args$method_name == "scrapper") {
  # scrapper --------------------------------------------------------------
  source(file.path(run_dir, "scrapper.R"))
  output_data <- run_scrapper(sce, n_cluster, args$n_comp, args$n_neig, 
                              args$n_hvg, args$filter, time, clustering_info)
} else if (args$method_name == "bixverse") {
  # scrapper --------------------------------------------------------------
  source(file.path(run_dir, "bixverse.R"))
  output_data <- run_bixverse(sce, n_cluster, args$n_comp, args$n_neig, 
                              args$n_hvg, args$filter, time, clustering_info)
}

# write outputs 
# timings -----------------------------------------------------------------
output_data$time <- lapply(output_data$time, as.numeric, units = "secs")
write_json(output_data$time, output[["methods_timings"]], 
           auto_unbox = TRUE, pretty = TRUE)
# clustering_info ---------------------------------------------------------
write_json(output_data$clustering_info, output[["methods_clustering_info"]],
           auto_unbox = TRUE, pretty = TRUE)
# clusters ----------------------------------------------------------------
write.table(data.frame(cell_id = output_data$cell_ids,
                       louvain = output_data$louvain, 
                       leiden = output_data$leiden),
            output[["methods_clusters"]], 
            sep = "\t", quote = FALSE, row.names = FALSE)
# PCAs --------------------------------------------------------------------
output_data$pca <- data.frame(cell_id = rownames(output_data$pca), 
                              output_data$pca)
colnames(output_data$pca) <- c("cell_id", paste0("PC", seq_len(ncol(output_data$pca) - 1)))
write.table(output_data$pca, output[["methods_pca"]], sep = "\t", quote = FALSE, row.names = FALSE)

cat(str(output_data))
cat(str(output_data$hvgs))
# HVGs --------------------------------------------------------------------
writeLines(output_data$hvgs, file(output[["methods_hvgs"]]))

