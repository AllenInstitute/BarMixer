#' Test if a hash distribution is unimodal
#'
#' This function uses the "SI" method provided by multimode::modetest() with B = 20.
#'
#' @param x a numeric vector of hash count values
#' @param p.cutoff a nuemric value indicating the p-value cutoff for a test of multimodality. If p.value > p.cutoff, the distribution is unimodal.
#'
#' @return a logical value. TRUE if x appears unimodal based on the p-value cutoff.
#'
is_unimodal <- function(x, p.cutoff = 0.2) {
  x <- log10(x + 1)
  p <- multimode::modetest(x, method = "SI", B = 20)$p.value
  p > p.cutoff
}

#' Select a cutoff for hash counts using automatic, k-means based splitting
#'
#' @param x a numeric vector of hash count values
#' @param use_median_cut a logical value indicating whether or not to use the median value of x to set a minimum threshold. Default is FALSE.
#' @param min_cut a numeric value for the minimum number of counts to consider. Default is 10.
#' @param min_fc a numeric value for the minimum expected fold change between the centers of the two clusters (positive and negative). Default is 2.
#' @param seed a value to use as a random seed for k-means clustering. Default is 3030.
#'
#' @return a numeric value to use as a cutoff. If no values are above min_cut, returns max(x), otherwise returns max(x) for non-passing values.
#' @export
#'
select_hash_cutoff <- function(x,
                               use_median_cut = FALSE,
                               min_cut = 10,
                               min_fc = 2,
                               seed = 3030) {

  assertthat::assert_that(class(x) == "numeric")
  assertthat::assert_that(class(use_median_cut) == "logical")
  assertthat::assert_that(length(use_median_cut) == 1)
  assertthat::is.number(min_cut)
  assertthat::assert_that(length(min_cut) == 1)
  assertthat::is.number(min_fc)
  assertthat::assert_that(length(min_fc) == 1)
  assertthat::is.number(seed)
  assertthat::assert_that(length(seed) == 1)

  # if no counts for this hash, return a cutoff of 0
  if(sum(x == 0) == length(x)) {
    return(0)
  }

  res <- rep(0, length(x))

  # if use_median_cut, use the median value to set the min cutoff
  # useful for cases when there are many low values
  # (usually the case for hashing)
  if(use_median_cut) {
    min_cut <- median(x[x > 0])
    if(is.na(min_cut)) {
      min_cut <- 0
    }
  }

  x_gt_cut <- x[x > min_cut]

  # If no values are greater than the cutoff, return max(x)
  if(length(x_gt_cut) == 0) {
    return(max(x))
  }

  # If x appears to be unimodal, set cutoff at mean - 2*sd
  if(is_unimodal(x)) {
    x_gt_cut <- log10(x_gt_cut + 1)
    log_cut <- mean(x_gt_cut) - 2*sd(x_gt_cut)
    res <- floor(10^log_cut)

    return(res)

  } else {

    if(length(x_gt_cut) > 2 & var(x_gt_cut) > 0) {
      set.seed(seed)
      km <- stats::kmeans(log10(x_gt_cut), centers = 2)
      cl <- km$cluster
      high_cl <- which(km$centers == max(km$centers))
      low_cl <- which(km$centers != max(km$centers))

      fc <- km$centers[high_cl] - km$centers[low_cl]
      if(fc > log10(min_fc)) {
        res[x > min_cut][cl == high_cl] <- 1
      }
    }

    if(sum(res) == 0) {
      max(x)
    } else {
      max(x[res == 0])
    }
  }

}

#' Binarize a vector of hash counts based on a cutoff value
#'
#' @param x a numeric vector of hash count values
#' @param cutoff a numeric value to use as a cutoff
#'
#' @return a binary vector with value 1 for all values in x > cutoff
#' @export
#'
binarize_hash <- function(x,
                          cutoff) {

  assertthat::assert_that(class(x) == "numeric")
  assertthat::is.number(cutoff)
  assertthat::assert_that(length(cutoff) == 1)

  res <- rep(0, length(x))

  res[x > cutoff] <- 1

  res
}

#' Add zero values for missing HTOs to a matrix
#'
#' @param mat A matrix or sparse matrix of hash tag oligo counts. Rows must be hashes, columns must be cells
#' @param valid_htos A character vector of valid hto sequences.
#'
#' @return a dgCMatrix of HTO counts with rows equal to valid_htos.
#' @export
#'
add_missing_hto_rows <- function(mat,
                                 valid_htos) {

  assertthat::assert_that(check_matrix(mat))
  assertthat::assert_that(class(valid_htos) == "character")

  if(class(mat) == "dgCMatrix") {
    mat <- as(mat, "matrix")
  }


  missing_htos <- setdiff(valid_htos, rownames(mat))
  if(length(missing_htos) > 0) {

    missing_mat <- matrix(0,
                          nrow = length(missing_htos),
                          ncol = ncol(mat))
    rownames(missing_mat) <- missing_htos
    colnames(missing_mat) <- colnames(mat)

    mat <- rbind(mat,
                 missing_mat)
  }

  mat <- mat[valid_htos,]

  as(mat,"dgCMatrix")
}

#' Binarize a matrix of hash counts
#'
#' If cutoff values are provided, they will be used. If not, cutoffs will be determined
#' using select_hash_cutoff()
#'
#' @param mat A matrix or sparse matrix of hash tag oligo counts. Rows must be hashes, columns must be cells
#' @param valid_htos (optional) a character vector of valid hto sequences. Default is NULL, which will use all sequences provided as rownames(mat).
#' @param use_median_cut a logical value indicating whether or not to use the median value of x to set a minimum threshold.
#' @param min_cut a numeric value for the minimum number of counts to consider for select_hash_cutoff(). Default is 10.
#' @param cutoff_vals (optional) a named numeric vector of cutoff values to apply. Length must be equal to nrow(mat), and names must exist in rownames(mat). Default is NULL, which will use select_hash_cutoff() to determine cutoff values.
#' @param expect_equal_loading a logical value indicating if equal loading is expected. This will try to correct overrepresented hashes only (not under-represented). Default is TRUE.
#' @param max_tries a numeric value indicating how many times to attempt adjustment if expect_equal_loading == TRUE to prevent runaway looping. Default is 5
#'
#' @return A list containing two objects:
#' \itemize{
#'   \item bmat: a binary sparse matrix (dgCMatrix) with dimensions equal to mat
#'   \item bsummary: a data.frame with summary values for each hash tag oligo: hash_barcode, cutoff, n_pos, n_neg, frac_pos, frac_neg
#' }
#' @export
#'
binarize_hash_matrix <- function(mat,
                                 valid_htos = NULL,
                                 use_median_cut = FALSE,
                                 min_cut = 10,
                                 cutoff_vals = NULL,
                                 expect_equal_loading = TRUE,
                                 max_tries = 5) {

  assertthat::assert_that(check_matrix(mat))
  assertthat::assert_that(class(expect_equal_loading) == "logical")
  assertthat::assert_that(class(use_median_cut) == "logical")
  assertthat::assert_that(length(use_median_cut) == 1)

  if("dgCMatrix" %in% class(mat)) {
    mat <- as(mat, "matrix")
  }

  if(!is.null(valid_htos)) {
    assertthat::assert_that(is.character(valid_htos))
    mat <- mat[rownames(mat) %in% valid_htos,]
  } else {
    valid_htos <- rownames(mat)
  }

  if(!is.null(cutoff_vals)) {
    assertthat::assert_that(is.numeric(cutoff_vals))
    assertthat::assert_that(names(cutoff_vals) %in% rownames(mat))
    assertthat::assert_that(length(cutoff_vals) == nrow(mat))

    cutoff_vals <- cutoff_vals[rownames(mat)]
  } else {
    assertthat::assert_that(class(min_cut) == "numeric")
    assertthat::assert_that(length(min_cut) == 1)

    cutoff_vals <- apply(mat,
                         1,
                         select_hash_cutoff,
                         use_median_cut = use_median_cut,
                         min_cut = min_cut)
  }

  bmat <- matrix(as.numeric(mat > cutoff_vals),
                 ncol = ncol(mat),
                 nrow = nrow(mat))
  rownames(bmat) <- rownames(mat)
  colnames(bmat) <- colnames(mat)

  if(expect_equal_loading) {
    # Try to correct cutoff if one appears overloaded
    # This is likely to occur if one hash has a high background population.
    n_pos <- rowSums(bmat)

    tries <- 0

    while(max(n_pos) > median(n_pos[n_pos > 0]) * 4 & tries < max_tries) {
      too_high <- which(n_pos == max(n_pos))
      cutoff_vals[too_high] <- select_hash_cutoff(mat[too_high,],
                                                  min_cut = cutoff_vals[too_high])

      bmat <- matrix(as.numeric(mat > cutoff_vals),
                     ncol = ncol(mat),
                     nrow = nrow(mat))
      rownames(bmat) <- rownames(mat)
      colnames(bmat) <- colnames(mat)

      n_pos <- rowSums(bmat)

      tries <- tries + 1
    }
  }

  bsummary <- data.frame(hto_barcode = rownames(bmat),
                         cutoff = cutoff_vals,
                         n_pos = rowSums(bmat),
                         n_neg = ncol(bmat) - rowSums(bmat),
                         n_below_threshold = rowSums(mat < min_cut),
                         frac_pos = rowSums(bmat) / ncol(bmat),
                         frac_neg = (ncol(bmat) - rowSums(bmat)) / ncol(bmat),
                         frac_below_threshold = rowSums(mat < min_cut) / ncol(bmat))

  missing_htos <- setdiff(valid_htos, rownames(bmat))
  if(length(missing_htos) > 0) {
    missing_summary <- data.frame(hto_barcode = missing_htos,
                                  cutoff = 0,
                                  n_pos = 0,
                                  n_neg = ncol(bmat),
                                  n_below_threshold = ncol(bmat),
                                  frac_pos = 0,
                                  frac_neg = 0,
                                  frac_below_threshold = 0)
    bsummary <- rbind(bsummary,
                      missing_summary)

  }

  bmat <- add_missing_hto_rows(bmat,
                               valid_htos)

  bsummary <- bsummary[match(valid_htos, bsummary$hto_barcode),]
  bsummary$frac_pos <- round(bsummary$frac_pos, 4)
  bsummary$frac_neg <- round(bsummary$frac_neg, 4)
  bsummary$frac_below_threshold <- round(bsummary$frac_below_threshold, 4)

  bmat <- bmat[valid_htos,]

  list(bmat = as(bmat, "dgCMatrix"),
       bsummary = bsummary)

}

#' Convert binary hash matrix results to categorical labels for each cell barcode
#'
#' @param bmat a matrix or sparse matrix with binarized hash scores, as generated by `binarize_hash_matrix()`. Rows must be hash barcodes, and columns must be cell barcodes.
#'
#' @return A list containing two objects:
#' \itemize{
#'   \item hash_category_table: a data.frame of categorical hash results: cell_barcode, hash_category, hash_sequence
#'   \item hash_summary: a data.frame with summary statistics for hashes in this run: hash_category, n_category, frac_category
#' }
#'
#' @export
#'
categorize_binary_hash_matrix <- function(bmat) {

  assertthat::assert_that(check_matrix(bmat))

  if(class(bmat) == "dgCMatrix") {
    bmat <- as(bmat, "matrix")
  }

  results <- apply(bmat, 2,
                   function(binary_hash_scores) {
                     hash_sum <- sum(binary_hash_scores)
                     if(hash_sum == 0) {
                       category <- "no_hash"
                       sequence <- NA
                     } else if(hash_sum == 1) {
                       category <- "singlet"
                       sequence <- rownames(bmat)[binary_hash_scores == 1]
                     } else if(hash_sum == 2) {
                       category <- "doublet"
                       sequence <- paste(rownames(bmat)[binary_hash_scores == 1],
                                         collapse = ";")
                     } else if(hash_sum > 2) {
                       category <- "multiplet"
                       sequence <- paste(rownames(bmat)[binary_hash_scores == 1],
                                         collapse = ";")
                     }
                     data.frame(hto_category = category,
                                hto_barcode = sequence)
                   })

  hash_category_table <- do.call("rbind",
                                 results)

  hash_category_table <- cbind(data.frame(cell_barcode = colnames(bmat)),
                               hash_category_table)

  category_count_table <- table(hash_category_table$hto_category)
  hash_summary <- data.frame(hto_category = names(category_count_table),
                             n_category = as.numeric(category_count_table),
                             frac_category = as.numeric(category_count_table) / sum(category_count_table))

  missing_summary <- data.frame(hto_category = "missing",
                                n_category = sum(is.na(hash_category_table$hash_category)),
                                frac_category = sum(is.na(hash_category_table$hash_category)) / nrow(hash_category_table))

  hash_summary <- rbind(hash_summary, missing_summary)

  required_categories <- c("no_hash","singlet","doublet","multiplet","missing")
  missing_categories <- setdiff(required_categories, hash_summary$hto_category)
  if(length(missing_categories) > 0) {
    zero_summary <- data.frame(hto_category = missing_categories,
                               n_category = 0,
                               frac_category = 0)
    hash_summary <- rbind(hash_summary,
                          zero_summary)
  }

  hash_summary <- hash_summary[match(required_categories, hash_summary$hto_category),]

  hash_summary$frac_category <- round(hash_summary$frac_category, 4)

  list(hash_category_table = hash_category_table,
       hash_summary = hash_summary)
}


#' Generate a summary table with the number and fraction of singlets with each HTO barcode
#'
#' @param hash_category_table A hash category table data.frame as produced by categorize_binary_hash_matrix()
#' @param valid_htos (optional) a character vector of valid hto sequences. Default is NULL, which will use unique barcodes in hash_category_table$hto_barcode.
#'
#' @return a data.frame with hto_barcode, n_singlets, and frac_singlets.
#' @export
#'
make_singlet_summary <- function(hash_category_table,
                                 valid_htos = NULL) {

  assertthat::assert_that(class(hash_category_table) == "data.frame")

  if(!is.null(valid_htos)) {
    assertthat::assert_that(is.character(valid_htos))
  } else {
    valid_htos <- unique(hash_category_table$hto_barcodes)
  }

  hash_category_table <- as.data.table(hash_category_table)
  hash_category_table <- hash_category_table[hash_category_table$hto_category == "singlet",]

  singlet_counts <- hash_category_table[,
                                        .(n_singlets = nrow(.SD)),
                                        by = "hto_barcode"]

  missing_htos <- setdiff(valid_htos, singlet_counts$hto_barcode)
  if(length(missing_htos) > 0) {
    missing_counts <- data.table(hto_barcode = missing_htos,
                                 n_singlets = 0)
    singlet_counts <- rbind(singlet_counts,
                            missing_counts)
  }

  singlet_counts$frac_singlets <- round(singlet_counts$n_singlets / sum(singlet_counts$n_singlets), 4)

  singlet_counts <- singlet_counts[match(valid_htos, singlet_counts$hto_barcode),]
  singlet_counts <- as.data.frame(singlet_counts)

  singlet_counts
}
