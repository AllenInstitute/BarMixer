#' Decide how many bits to use to store integer values in an HDF5 Dataset object.
#'
#' @param x a numeric vector
#'
#' @return an integer value. One of 16L, 32L, 64L, or NULL if larger than 64-bit.
#' @export
#'
choose_integer_bits <- function(x) {
  max_val <- max(x, na.rm = TRUE)
  
  if(max_val < 2^16) {
    16L
  } else if(max_val < 2^32) {
    32L
  } else if(max_val < 2^64) {
    64L
  } else {
    NULL
  }
  
}

#' Select a reasonable chunk size for an HDF5 dataset object
#'
#' If x has length > 1e6, chunks will be 1e5.
#' If x has length > 1e2, chunks will be one log10 smaller than the length of x.
#' If x has length <= 1e2, chunk size is the length of x.
#'
#' @param x A vector to store in an HDF5 file
#'
#' @return a numeric value containing a suggested chunk size
#' @export
#'
choose_chunk_size <- function(x) {
  x_len <- length(x)
  x_logs <- floor(log10(x_len))
  
  if(x_logs > 6) {
    10^5
  } else if (x_logs > 2) {
    10^(x_logs - 1)
  } else {
    length(x)
  }
}

#' Generate a set of attributes based on 10x Genomics defaults
#'
#' @param library_ids Library IDs to store as attributes. if NULL (default), will generate a random ID using date and ids::proquint().
#'
#' @return a list of attributes
#' @export
#'
h5_attr_list <- function(library_ids = NULL) {
  
  attr_list <- list(chemistry_description = "Single Cell 3' v3",
                    filetype = "matrix",
                    library_ids = paste(Sys.Date(),
                                        ids::proquint(1),
                                        sep = "-"),
                    original_gem_groups = 1,
                    version = 2)
  if(!is.null(library_ids)) {
    assertthat::assert_that(class(library_ids) == "character")
    attr_list$library_ids <- library_ids
  }
  
  attr_list
}

#' Write an h5_list, as created by rhdf5::h5dump(), to an .h5 file
#'
#' @param h5_list a list object, e.g. a list created by rhdf5::h5dump()
#' @param h5_file a character object specifying the location of a .h5 file to write to.
#' @param overwrite a logical value specifying whether or not to overwrite an existing .h5 file. Default is FALSE.
#' @param h5_handle an existing h5_handle created by H5Fopen(). Used for recursion. The default (NULL) should usually be used.
#' @param h5_target a base location within the HDF5 file to write to. Mainly used for recursion. The default ("/") should usually be used.
#' @param h5_attributes a list of attributes to add to an .h5 file to try to imitate 10x Genomics outputs. If NULL (default), will be skipped. "tenx" uses in-built data from 'h5_attr_list()'.
#' @param library_ids a character vector of library ids to add to attributes. Only used if h5_attributes != NULL. Default is NULL.
#'
#' @return Writes a file; no return to R.
#' @export
#'
write_h5_list <- function(h5_list,
                          h5_file,
                          overwrite = FALSE,
                          addon = FALSE,
                          h5_handle = NULL,
                          h5_target = "/",
                          h5_attributes = NULL,
                          library_ids = NULL) {
  
  assertthat::assert_that(is.list(h5_list))
  assertthat::assert_that(is.character(h5_file))
  assertthat::assert_that(length(h5_file) == 1)
  
  if(!is.null(h5_attributes)) {
    if(h5_attributes  == "tenx") {
      h5_attributes <- BarMixer::h5_attr_list()
    }
  }
  
  if(!is.null(library_ids)) {
    h5_attributes$library_ids <- library_ids
  }
  
  # Make sure the HDF5 file connection is closed if the function
  # exits due to an error.
  on.exit(expr = {
    if(h5_target == "/") {
      rhdf5::h5closeAll()
    }
  })
  
  if(is.null(h5_handle)) {
    if(!overwrite & !addon) {
      # If we aren't overwriting or adding, check for file and halt if it exists.
      # If it doesn't exist, create it.
      
      if(file.exists(h5_file)) {
        stop(paste(h5_file, "already exists."))
      } else {
        rhdf5::H5Fcreate(h5_file)
      }
      
    } else if(overwrite) {
      # If we're overwriting, remove the old file and make a new one.
      # Partial changes are performed by addon
      
      # If overwrite and addon, overwrite takes precedence.
      
      if(file.exists(h5_file)) {
        file.remove(h5_file)
      }
      
      rhdf5::H5Fcreate(h5_file)
      
    } else if(addon) {
      # If addon, we don't need to create or remove the file.
      # Keeping this condition to help with reasoning even though
      # nothing happens here.
      # We'll open below.
    } else {
      rhdf5::H5Fcreate(h5_file)
    }
    
    h5_handle <- rhdf5::H5Fopen(h5_file)
  }
  
  # Add file attributes to match cellranger
  if(!is.null(h5_attributes) & !addon) {
    if(h5_target != "/" & class(h5_attributes) == "list") {
      base_obj <- rhdf5::H5Dopen(h5_handle, "/")
      
      for(i in 1:length(h5_attributes)) {
        rhdf5::h5writeAttribute(h5_attributes[[i]],
                                base_obj,
                                names(h5_attributes)[i])
      }
      
      rhdf5::H5Dclose(base_obj)
    }
    h5_attributes <- NULL
  }
  
  
  h5_names <- names(h5_list)
  if(length(h5_names) > 0) {
    for(h5_name in h5_names) {
      new_object <- paste0(h5_target, h5_name)
      
      # Correct 1d arrays to vectors for storage
      if(class(h5_list[[h5_name]]) == "array") {
        if(length(dim(h5_list[[h5_name]])) == 1) {
          h5_list[[h5_name]] <- as.vector(h5_list[[h5_name]])
        }
      }
      
      # Correct numeric arrays to integers if they're all whole numbers
      if(class(h5_list[[h5_name]]) == "numeric") {
        
        x <- as.integer(h5_list[[h5_name]])
        if(isTRUE(all.equal(h5_list[[h5_name]], x, check.attributes = FALSE))) {
          h5_list[[h5_name]] <- x
        }
      }
      
      if(class(h5_list[[h5_name]]) == "list") {
        rhdf5::h5createGroup(h5_handle,
                             group = new_object)
        # Recurse function to write children of list
        BarMixer::write_h5_list(h5_list[[h5_name]],
                                h5_file = h5_file,
                                h5_handle = h5_handle,
                                h5_target = paste0(new_object,"/"),
                                h5_attributes = h5_attributes,
                                addon = addon)
      } else if(class(h5_list[[h5_name]]) == "numeric") {
        rhdf5::h5createDataset(h5_handle,
                               dataset = new_object,
                               dims = list(length(h5_list[[h5_name]])),
                               chunk = BarMixer::choose_chunk_size(h5_list[[h5_name]]),
                               storage.mode = storage.mode(h5_list[[h5_name]]))
        
        rhdf5::h5write(obj = h5_list[[h5_name]],
                       file = h5_handle,
                       name = new_object)
        
      } else if(class(h5_list[[h5_name]]) == "integer") {
        bits <- BarMixer::choose_integer_bits(h5_list[[h5_name]])
        h5_type <- paste0("H5T_NATIVE_UINT",bits)
        
        rhdf5::h5createDataset(h5_handle,
                               dataset = new_object,
                               dims = list(length(h5_list[[h5_name]])),
                               chunk = BarMixer::choose_chunk_size(h5_list[[h5_name]]),
                               H5type = h5_type)
        
        rhdf5::h5write(obj = h5_list[[h5_name]],
                       file = h5_handle,
                       name = new_object)
        
      } else if(class(h5_list[[h5_name]]) == "character") {
        h5_list[[h5_name]][is.na(h5_list[[h5_name]])] <- "NA"
        
        rhdf5::h5createDataset(h5_handle,
                               dataset = new_object,
                               dims = list(length(h5_list[[h5_name]])),
                               chunk = BarMixer::choose_chunk_size(h5_list[[h5_name]]),
                               storage.mode = storage.mode(h5_list[[h5_name]]),
                               size = max(nchar(h5_list[[h5_name]])) + 1)
        
        rhdf5::h5write(obj = h5_list[[h5_name]],
                       file = h5_handle,
                       name = new_object)
      }
    }
  }
  
}

#' Add a transposed version of a matrix to an .h5 file
#'
#' Useful for applications that need feature-indexed, rather than observation-indexed values.
#'
#' @param h5_file An existing .h5 file.
#' @param in_target The name of the target matrix. Default is "matrix".
#' @param in_feature_names The feature names to use in the transposed matrix. Default is "id".
#' @param in_sample_names The observation names to use in the transposed matrix. Default is "barcodes".
#' @param out_target The name of the output matrix to store. Default is "transposed_" prefixed to the in_target param (e.g. "transposed_matrix").
#'
#' @return No return
#' @export
#'
add_h5_transposed_matrix <- function(h5_file,
                                     in_target = "matrix",
                                     in_feature_names = "id",
                                     in_sample_names = "barcodes",
                                     out_target = paste0("transposed_",target)) {
  
  assertthat::assert_that(is.character(h5_file))
  assertthat::assert_that(length(h5_file) == 1)
  assertthat::assert_that(file.exists(h5_file))
  
  mat <- BarMixer::read_h5_dgCMatrix(h5_file,
                                     target = target,
                                     feature_names = in_feature_names,
                                     sample_names = in_sample_names)
  
  mat <- list(Matrix::t(mat))
  names(mat) <- paste0(out_target,"_dgCMatrix")
  
  mat <- BarMixer::h5_list_convert_from_dgCMatrix(mat,
                                                  target = out_target)
  
  BarMixer::write_h5_list(mat,
                          h5_file = h5_file,
                          overwrite = FALSE,
                          addon = TRUE)
}
