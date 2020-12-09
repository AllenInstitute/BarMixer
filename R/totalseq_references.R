#' Load BioLegend TotalSeq A reference table for human hash tag oligos (HTO)
#'
#' Note: there is no Hashtag 11 in the human HTO set.
#'
#' @return a data.table with reference information
#' @export
#'
#' @importFrom data.table fread
#'
totalseq_a_human <- function() {
  dt <- fread(
    system.file("reference/TotalSeqA_human_barcodes.csv",
                package = "HTOparser")
    )
  dt[order(dt$hto_id),]
}

#' Load BioLegend TotalSeq A reference table for mouse hash tag oligos (HTO)
#'
#' @return a data.table with reference information
#' @export
#'
#' @importFrom data.table fread
#'
totalseq_a_mouse <- function() {
  dt <- fread(
    system.file("reference/TotalSeqA_mouse_barcodes.csv",
                package = "HTOparser")
  )
  dt[order(dt$hto_id),]
}

#' Load BioLegend TotalSeq A reference table for biotin hash tag oligos (HTO)
#'
#' @return a data.table with reference information
#' @export
#'
#' @importFrom data.table fread
#'
totalseq_a_biotin <- function() {
  dt <- fread(
    system.file("reference/TotalSeqA_biotin_barcodes.csv",
                package = "HTOparser")
  )
  dt[order(dt$hto_id),]
}
