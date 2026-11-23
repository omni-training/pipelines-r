suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(bixverse)
})

run_bixverse <- function(
  sce, n_cluster, n_comp = 50, n_neig = 15, n_hvg = 1000,
  filter = c("manual", "auto"), time, clustering_info
) {
  filter <- match.arg(filter)
  nthreads <- 1
  # bixverse runs its parallel work in Rust; this caps the thread pool if it
  # has not already been started in this R session
  Sys.setenv(RAYON_NUM_THREADS = nthreads)

  # read the manual QC thresholds before the data leaves the SCE
  qc_thresholds <- metadata(sce)$qc_thresholds

  # bixverse takes cell and gene IDs from colData / rowData
  colData(sce)$cell_id <- colnames(sce)
  rowData(sce)$gene_id <- rownames(sce)

  dir_data <- tempfile("bixverse_")
  dir.create(dir_data, recursive = TRUE)
  on.exit(unlink(dir_data, recursive = TRUE, force = TRUE), add = TRUE)

  # load data + normalization ####
  # bixverse normalises the counts while loading them into its on-disk store,
  # so this step stands in for scrapper's normalization step
  start_time <- Sys.time()
  sc_object <- load_sce(
    object = SingleCells(dir_data = dir_data),
    sce = sce,
    assay_name = "counts",
    # minimal thresholds so the loader keeps (almost) all cells and genes;
    # the real filtering happens below
    sc_qc_param = params_sc_min_quality(
      min_unique_genes = 1L,
      min_lib_size = 1L,
      min_cells = 1L
    ),
    .verbose = FALSE
  )
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Normalization (load + normalise). Time Elapsed:", time_elapsed))
  time$normalization <- time_elapsed

  #### 1. find mithocondial genes ####
  start_time <- Sys.time()
  var <- get_sc_var(sc_object)
  gs_of_interest <- list(
    MT = var[grepl("^MT-", gene_id), gene_id]
  )
  sc_object <- gene_set_proportions_sc(
    sc_object,
    gs_of_interest,
    streaming = FALSE,
    .verbose = FALSE
  )
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Find mithocondrial genes. Time Elapsed:", time_elapsed))
  time$find_mit_gene <- time_elapsed

  # filter data ####
  start_time <- Sys.time()
  qc_df <- sc_object[[c("cell_id", "lib_size", "nnz", "MT")]]
  if (filter == "manual") {
    qc <- qc_thresholds
    mt_percent <- qc_df$MT * 100
    write(paste0("mt pcts: ", paste0(round(head(mt_percent, 3), 8), collapse = ";")), stderr())
    write(paste0("cells detected: ", paste0(round(head(qc_df$nnz, 3), 8), collapse = ";")), stderr())
    write(paste0("cells sum: ", paste0(round(head(qc_df$lib_size, 3), 8), collapse = ";")), stderr())
    keep <- qc_df$nnz >= qc[qc$metric == "nFeature", "min"] &
      qc_df$nnz <= qc[qc$metric == "nFeature", "max"] &
      mt_percent < qc[qc$metric == "percent.mt", "max"] &
      qc_df$lib_size <= qc[qc$metric == "nCount", "max"]
  } else {
    # MAD-based outlier detection, as in the bixverse PBMC vignette
    metrics <- list(
      log10_lib_size = log10(qc_df$lib_size),
      log10_nnz = log10(qc_df$nnz),
      MT = qc_df$MT
    )
    directions <- c(
      log10_lib_size = "twosided",
      log10_nnz = "twosided",
      MT = "above"
    )
    cell_qc <- run_cell_qc(
      metrics = metrics,
      cells_to_keep = get_cells_to_keep(sc_object),
      directions = directions,
      threshold = 3
    )
    sc_object[["outlier"]] <- cell_qc$combined
    keep <- !cell_qc$combined
  }
  write(paste0("cells before: ", nrow(qc_df)), stderr())
  sc_object <- set_cells_to_keep(sc_object, qc_df[keep, cell_id])
  write(paste0("cells after: ", length(get_cells_to_keep(sc_object))), stderr())
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Filter data. Time Elapsed:", time_elapsed))
  time$filter <- time_elapsed

  # Identification of highly variable features (feature selection) ####
  start_time <- Sys.time()
  sc_object <- find_hvg_sc(
    object = sc_object,
    hvg_no = as.integer(n_hvg),
    .verbose = FALSE
  )
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Identification of highly variable features. Time Elapsed:", time_elapsed))
  time$hvg <- time_elapsed

  # PCA ####
  start_time <- Sys.time()
  sc_object <- calculate_pca_sc(
    object = sc_object,
    no_pcs = as.integer(n_comp),
    sparse_svd = TRUE
  )
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("PCA. Time Elapsed:", time_elapsed))
  time$pca <- time_elapsed

  # t-sne ####
  start_time <- Sys.time()
  sc_object <- tsne_sc(sc_object, .verbose = FALSE)
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("t-sne. Time Elapsed:", time_elapsed))
  time$t_sne <- time_elapsed

  # umap ####
  start_time <- Sys.time()
  # use_knn = FALSE: neighbours are computed from the PCA, as in runUmap()
  sc_object <- umap_sc(sc_object, use_knn = FALSE, seed = 42L, .verbose = FALSE)
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("UMAP. Time Elapsed:", time_elapsed))
  time$umap <- time_elapsed

  # louvain ####
  start_time <- Sys.time()
  sc_object <- find_neighbours_sc(
    object = sc_object,
    neighbours_params = params_sc_neighbours(
      knn = list(k = as.integer(n_neig))
    ),
    .verbose = FALSE
  )
  louvain_search <- binary_search(
    sc_object,
    do_clustering = function(spe, resolution) {
      spe <- find_clusters_sc(
        spe,
        cluster_algorithm = "louvain",
        res = resolution,
        name = "louvain",
        .verbose = FALSE
      )
      get_sc_obs(spe, filtered = TRUE)[["louvain"]]
    },
    extract_nclust = function(result) length(unique(result)),
    n_clust_target = n_cluster
  )
  louvain_clustering <- louvain_search$result
  clustering_info$resolutions$louvain <- louvain_search$resolution
  clustering_info$num_runs$louvain <- louvain_search$num_runs
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  avg_time_elapsed <- time_elapsed / louvain_search$num_runs
  print(paste("Louvain clusterings. Total Search Time Elapsed:", time_elapsed))
  print(paste("Louvain clusterings. Average Time per Run:", avg_time_elapsed))
  print(paste("Louvain resolution:", clustering_info$resolutions$louvain))
  print(paste("Louvain runs:", clustering_info$num_runs$louvain))
  time$louvain <- avg_time_elapsed

  # leiden ####
  start_time <- Sys.time()
  sc_object <- find_neighbours_sc(
    object = sc_object,
    neighbours_params = params_sc_neighbours(
      knn = list(k = as.integer(n_neig))
    ),
    .verbose = FALSE
  )
  leiden_search <- binary_search(
    sc_object,
    do_clustering = function(spe, resolution) {
      spe <- find_clusters_sc(
        spe,
        cluster_algorithm = "leiden",
        res = resolution,
        name = "leiden",
        .verbose = FALSE
      )
      get_sc_obs(spe, filtered = TRUE)[["leiden"]]
    },
    extract_nclust = function(result) length(unique(result)),
    n_clust_target = n_cluster
  )
  leiden_clustering <- leiden_search$result
  clustering_info$resolutions$leiden <- leiden_search$resolution
  clustering_info$num_runs$leiden <- leiden_search$num_runs
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  avg_time_elapsed <- time_elapsed / leiden_search$num_runs
  print(paste("Leiden clusterings. Total Search Time Elapsed:", time_elapsed))
  print(paste("Leiden clusterings. Average Time per Run:", avg_time_elapsed))
  print(paste("Leiden resolution:", clustering_info$resolutions$leiden))
  print(paste("Leiden runs:", clustering_info$num_runs$leiden))
  time$leiden <- avg_time_elapsed

  cell_ids <- get_sc_obs(sc_object, filtered = TRUE)[["cell_id"]]
  pca <- get_pca_factors(sc_object)
  rownames(pca) <- cell_ids

  return(list(
    pca = pca,
    hvgs = unname(get_gene_names_from_idx(sc_object, get_hvg(sc_object))),
    cell_ids = cell_ids,
    time = time,
    clustering_info = clustering_info,
    leiden = leiden_clustering,
    louvain = louvain_clustering
  ))
}


