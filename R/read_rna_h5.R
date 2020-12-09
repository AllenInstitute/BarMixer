#' Read the /matrix from a .h5 file as a sparse matrix
#'
#' @param h5_file the path to an .h5 file in 10x Genomics format
#' @param target a character object specifying the target matrix within the file. Default is "matrix".
#' @param feature_names a character object specifying whether to use "id" or "name" for row.names. Default is "id".
#' @param sample_names a character object specifying which values to use for col.names. If "barcodes", will use /target/barcodes. Other values will be read from /target/observations/
#' @param index1 a logical object specifying whether index vectors should start with 0 (FALSE) or 1 (TRUE). Default is TRUE.
#'
#' @return a dgCMatrix of gene expression values.
#' @export
#'
read_h5_dgCMatrix <- function(h5_file,
                              target = "matrix",
                              feature_names = "id",
                              sample_names = "barcodes",
                              index1 = TRUE) {

  assertthat::assert_that(is.character(h5_file))
  assertthat::assert_that(length(h5_file) == 1)

  assertthat::assert_that(is.character(target))
  assertthat::assert_that(length(target) == 1)

  if(grepl("^/",target)) {
    target <- sub("^/","",target)
  }

  assertthat::assert_that(is.character(feature_names))
  assertthat::assert_that(length(feature_names) == 1)

  assertthat::assert_that(is.character(sample_names))
  assertthat::assert_that(length(sample_names) == 1)

  if(!file.exists(h5_file)) {
    stop(paste(h5_file, "does not exist."))
  }

  # Make sure the HDF5 file connection is closed if the function
  # exits due to an error.
  on.exit(expr = {
    rhdf5::H5Fclose(h5_handle)
  })

  feature_names <- match.arg(arg = feature_names,
                             choices = c("id","name"))

  h5_handle <- rhdf5::H5Fopen(h5_file)

  if(sample_names == "barcodes") {
    colname_target <- paste0("/", target, "/barcodes")
  } else {
    colname_target <- paste0("/", target, "/observations/", sample_names)
  }

  if(index1) {
    mat <- Matrix::sparseMatrix(x = rhdf5::h5read(h5_handle, paste0("/",target,"/data")),
                                i = rhdf5::h5read(h5_handle, paste0("/",target,"/indices")) + 1,
                                p = rhdf5::h5read(h5_handle, paste0("/",target,"/indptr")),
                                index1 = index1,
                                dims = rhdf5::h5read(h5_handle, paste0("/",target,"/shape")),
                                dimnames = list(as.vector(rhdf5::h5read(h5_handle, paste0("/",target,"/features/",feature_names))),
                                                as.vector(rhdf5::h5read(h5_handle, colname_target))
                                )
    )
  } else {
    mat <- Matrix::sparseMatrix(x = rhdf5::h5read(h5_handle, paste0("/",target,"/data")),
                                i = rhdf5::h5read(h5_handle, paste0("/",target,"/indices")),
                                p = rhdf5::h5read(h5_handle, paste0("/",target,"/indptr")),
                                index1 = index1,
                                dims = rhdf5::h5read(h5_handle, paste0("/",target,"/shape")),
                                dimnames = list(as.vector(rhdf5::h5read(h5_handle, paste0("/",target,"/features/",feature_names))),
                                                as.vector(rhdf5::h5read(h5_handle, colname_target))
                                )
    )
  }

  mat
}

#' Read .h5 Cell Metadata
#'
#' @param h5_file the path to an .h5 file in 10x Genomics format
#' @param target A matrix object in the .h5 file with a /barcodes object and/or a /target/observations/ sub-group. Default is "matrix".
#'
#' @return A data.frame containing all feature metadata found in /target/barcodes and /target/observations/
#' @export
#'
read_h5_cell_meta <- function(h5_file,
                              target = "matrix") {

  assertthat::assert_that(is.character(h5_file))
  assertthat::assert_that(length(h5_file) == 1)

  target <- ifelse(grepl("^/",target),
                   target,
                   paste0("/",target))

  h5_contents <- H5weaver::h5ls(h5_file)
  target_contents <- h5_contents[grepl(paste0("^",target), h5_contents$group),]

  h5_meta_targets <- character()

  target_bcs <- paste0(target, "/barcodes")

  if(target_bcs %in% target_contents$full_name) {
    h5_meta_targets <- c(h5_meta_targets,
                         target_bcs)
  }

  target_obs <- paste0(target, "/observations")

  if(target_obs %in% target_contents$full_name) {
    h5_meta_targets <- c(h5_meta_targets,
                         target_contents$full_name[target_contents$group == target_obs])
  }

  if(length(h5_meta_targets) > 0) {
    meta_list <- lapply(h5_meta_targets,
                        function(h5_meta_target) {
                          rhdf5::h5read(h5_file,
                                        h5_meta_target)
                        })
    rhdf5::h5closeAll()

    names(meta_list) <- sub(".+/","",h5_meta_targets)

    meta_list <- strip_1d_array_recursive(meta_list)
    meta_list <- convert_char_na_recursive(meta_list)

    df <- as.data.frame(meta_list,
                        stringsAsFactors = FALSE)

    df
  } else {
    stop("No cell metadata found in h5_file.")
  }
}

#' Read .h5 Feature Metadata
#'
#' @param h5_file the path to an .h5 file in 10x Genomics format
#' @param target A matrix object in the .h5 file with a /features/ sub-group. Default is "matrix".
#'
#' @return a data.frame containing all feature metadata found in /target/features/
#' @export
read_h5_feature_meta <- function(h5_file,
                                 target = "matrix") {
  assertthat::assert_that(is.character(h5_file))
  assertthat::assert_that(length(h5_file) == 1)

  target <- ifelse(grepl("^/",target),
                   target,
                   paste0("/",target))

  h5_contents <- H5weaver::h5ls(h5_file)
  target_contents <- h5_contents[grepl(paste0("^",target), h5_contents$group),]

  h5_meta_targets <- character()

  target_feat <- paste0(target, "/features")

  if(target_feat %in% target_contents$full_name) {
    h5_meta_targets <- c(h5_meta_targets,
                         target_contents$full_name[target_contents$group == target_feat])
  }

  h5_meta_names <- sub(".+/","",h5_meta_targets)
  h5_meta_targets <- h5_meta_targets[!grepl("^_",h5_meta_names)]
  h5_meta_names <- h5_meta_names[!grepl("^_",h5_meta_names)]

  if(length(h5_meta_targets) > 0) {
    meta_list <- lapply(h5_meta_targets,
                        function(h5_meta_target) {
                          rhdf5::h5read(h5_file,
                                        h5_meta_target)
                        })
    rhdf5::h5closeAll()

    names(meta_list) <- h5_meta_names

    meta_list <- strip_1d_array_recursive(meta_list)
    meta_list <- convert_char_na_recursive(meta_list)

    df <- as.data.frame(meta_list,
                        stringsAsFactors = FALSE)

    df
  } else {
    stop("No cell metadata found in h5_file.")
  }
}

#' Read .h5 directly into a Seurat Object
#'
#' By default, matrix data will be stored in the "RNA" assay.
#'
#' CITE-seq data will be automatically detected and stored in the "ADT" assay.
#'
#' @param h5_file the path to an .h5 file in 10x Genomics format
#' @param target A matrix object in the .h5 file with a /features/ sub-group. Default is "matrix".
#' @param feature_names a character object specifying whether to use "id" or "name" for row.names. Default is "id".
#' @param ... Additional parameters passed to \code{\link{Seurat::createSeuratObject}}
#'
#' @return a Seurat Class object
#' @export
read_h5_seurat <- function(h5_file,
                           target = "matrix",
                           feature_names = "id",
                           ...) {

  if(!requireNamespace("Seurat", versionCheck = list(op = ">=", version = "3.1.0"))) {
    stop("Can't find the Seurat package. Please install with install.packages(\"Seurat\")")
  }

  assertthat::assert_that(is.character(h5_file))
  assertthat::assert_that(length(h5_file) == 1)
  assertthat::assert_that(typeof(feature_names) == "character")
  assertthat::assert_that(length(feature_names) == 1)
  assertthat::assert_that(feature_names %in% c("id","name"))

  mat <- read_h5_dgCMatrix(h5_file,
                           target = target,
                           feature_names = feature_names)

  cell_meta <- read_h5_cell_meta(h5_file,
                                 target = target)
  rownames(cell_meta) <- cell_meta$barcodes

  rownames(cell_meta) <- cell_meta$barcodes

  feat_meta <- read_h5_feature_meta(h5_file,
                                    target = target)

  cite <- FALSE

  # Check for CITE-seq data
  if("feature_type" %in% names(feat_meta)) {
    if("Antibody Capture" %in% feat_meta$feature_type) {
      cite <- TRUE
    }
  }

  if(cite) {
    cite_feat <- feat_meta[feat_meta$feature_type == "Antibody Capture",]
    feat_meta <- feat_meta[feat_meta$feature_type != "Antibody Capture",]

    cite_mat <- mat[cite_feat$id,]
    mat <- mat[feat_meta$id,]
  }

  so <- Seurat::CreateSeuratObject(counts = mat,
                                   meta.data = cell_meta,
                                   ...)
  if(cite) {
    so[["ADT"]] <- Seurat::CreateAssayObject(counts = cite_mat)
  }

  so
}

#' Read .h5 directly into a SingleCellExperiment object
#'
#' By default, matrix data will be stored in the "counts" assay.
#'
#' CITE-seq data will be automatically detected and stored in the altExp called "ADT"
#'
#' @param h5_file the path to an .h5 file in 10x Genomics format
#' @param target A matrix object in the .h5 file with a /features/ sub-group. Default is "matrix".
#' @param ... Additional parameters passed to \code{\link{SingleCellExperiment::SingleCellExperiment()}}
#'
#' @return a SingleCellExperiment Class object
#' @export
read_h5_sce <- function(h5_file,
                        target = "matrix",
                        ...) {

  if(!requireNamespace("SingleCellExperiment", versionCheck = list(op = ">=", version = "1.8.0"))) {
    stop("Can't find the SingleCellExperiment package. Please install with BiocManater::install(\"SingleCellExperiment\")")
  }

  assertthat::assert_that(is.character(h5_file))
  assertthat::assert_that(length(h5_file) == 1)

  mat <- read_h5_dgCMatrix(h5_file,
                           target = target)

  cell_meta <- read_h5_cell_meta(h5_file,
                                 target = target)

  feat_meta <- read_h5_feature_meta(h5_file,
                                    target = target)

  cite <- FALSE

  # Check for CITE-seq data
  if("feature_type" %in% names(feat_meta)) {
    if("Antibody Capture" %in% feat_meta$feature_type) {
      cite <- TRUE
    }
  }

  if(cite) {
    cite_feat <- feat_meta[feat_meta$feature_type == "Antibody Capture",]
    feat_meta <- feat_meta[feat_meta$feature_type != "Antibody Capture",]

    cite_mat <- mat[cite_feat$id,]
    mat <- mat[feat_meta$id,]
  }

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = mat),
    colData = cell_meta,
    rowData = feat_meta,
    ...
  )

  if(cite) {
    cite_se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = cite_mat)
    )
    SingleCellExperiment::altExp(sce, "ADT") <- cite_se
  }

  sce
}

#' List objects in an HDF5 file
#'
#' This is a wrapper around rhdf5::h5ls() that adds the full name of each object.
#'
#' @param ... parameters passed to rhdf5::h5ls()
#'
#' @return a data.frame listing the contents of the target hdf5 file.
#' @export
#'
h5ls <- function(...) {
  df <- rhdf5::h5ls(...)

  df$full_name <- paste0(df$group, "/" ,df$name)
  df$full_name <- sub("//","/",df$full_name)

  df <- df[,c("full_name","group","name","otype","dclass","dim")]

  df
}

#' Check if an object exists in an HDF5 file
#'
#' @param h5_file A character object specifying the path to a .h5 file
#' @param target A character object specifying the name of the object to check (e.g. "/matix/features/id")
#'
#' @return a logical value
#' @export
h5exists <- function(h5_file,
                     target) {

  assertthat::assert_that(is.character(h5_file))
  assertthat::assert_that(length(h5_file) == 1)
  assertthat::assert_that(file.exists(h5_file))

  assertthat::assert_that(is.character(target))
  assertthat::assert_that(length(target) == 1)

  df <- H5weaver::h5ls(h5_file)

  target %in% df$full_name
}

#' Convert all 1D Arrays in a list object to vectors recursively
#'
#' @param x a list object
#'
#' @return a list object with all 1D arrays converted.
#' @export
#'
strip_1d_array_recursive <- function(x) {
  assertthat::assert_that(class(x) == "list")

  if(length(x) > 0) {
    for(n in seq_along(x)) {
      if(class(x[[n]]) == "list") {
        x[[n]] <- strip_1d_array_recursive(x[[n]])
      } else if(class(x[[n]]) == "array" & length(dim(x[[n]] == 1))) {
        x[[n]] <- as.vector(x[[n]])
      }
    }
  }

  x
}


#' Convert all "NA" character values to actual NAs recursively
#'
#' @param x a list object
#'
#' @return a list object with "NA"s converted to NA
#' @export
#'
convert_char_na_recursive <- function(x) {
  assertthat::assert_that(class(x) == "list")

  if(length(x) > 0) {
    for(n in seq_along(x)) {
      if(class(x[[n]]) == "list") {
        x[[n]] <- convert_char_na_recursive(x[[n]])
      } else if(class(x[[n]]) == "character") {
        x[[n]] <- convert_char_na(x[[n]])
      }
    }
  }

  x
}

#' Dump all objects from an HDF5 file to a list.
#'
#' This is a wrapper around rhdf5::h5dump() that converts all 1D arrays to vectors
#' and correctly handles NA values.
#'
#' @param ... parameters passed to rhdf5::h5dump()
#'
#' @return a list object with the contents of the target HDF5 file
#' @export
#'
h5dump <- function(...) {

  h5_list <- rhdf5::h5dump(...)
  h5_list <- strip_1d_array_recursive(h5_list)
  h5_list <- convert_char_na_recursive(h5_list)

  h5_list
}

#' Get dimensions of an object in an HDF5 file
#'
#' @param h5_file The path to an HDF5 file
#' @param name The full name of an object in the HDF5 file.
#'
#' @return a numeric vector containing the dimensions of the target object.
#' @export
#'
h5dims <- function(h5_file,
                   name) {

  h5_contents <- H5weaver::h5ls(h5_file)

  d <- h5_contents$dim[h5_contents$full_name == name]
  d <- unlist(strsplit(d, split = ","))

  as.numeric(d)
}
