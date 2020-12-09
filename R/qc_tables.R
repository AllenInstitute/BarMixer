#' Render a nicely-formatted table for QC reports
#'
#' @param df A data.frame to format
#'
#' @return an HTML/JavaScript DT object.
#' @export
#'
qc_table <- function(df) {

  assertthat::assert_that(sum(class(df) %in% c("data.frame","data.table")) > 0)

  if("data.table" %in% class(df)) {
    df <- as.data.frame(df)
  }

  # DOM for table
  # Based on bootstrap. Each row has 12 column units
  # l length input control (Show N Rows)
  # f filtering input (Search box)
  # r processing elements
  # t table
  # i info summary (which rows are displayed)
  # B buttons
  # p pagination
  row1 <- "<'row'<'col-md-6'l><'col-md-6'fr>>"
  row2 <- "<'row'<'col-md-12't>>"
  row3 <- "<'row'<'col-md-12'i>>"
  row4 <- "<'row'<'col-md-6'B><'col-md-6'p>>"

  dom <- paste0(row1, row2, row3, row4)

  DT::datatable(df,
                class = "compact table-striped",
                extensions = "Buttons",
                style = "bootstrap",
                options = list(dom = dom,
                               buttons = c('copy', 'csv', 'excel'),
                               lengthMenu = c(15, 25, 50),
                               scrollX = TRUE),
                rownames = FALSE)
}


#' Flip the orientation of a data.frame and return as a list.
#'
#' Mainly for use in formatting tables for storage in JSON
#'
#' @param df a data.frame or data.frame-like object
#' @param name_col (optional) a column in the data.frame used to name each row in the list values
#'
#' @return a nested list object.
#' @export
#'
qc_flip <- function(df,
                    name_col = NULL) {

  if(!is.null(name_col)) {
    apply_cols <- setdiff(1:ncol(df), which(names(df) == name_col))
  } else {
    apply_cols <- 1:ncol(df)
  }

  out_list <- lapply(1:nrow(df),
                     function(row_id) {
                       values <- lapply(apply_cols,
                                        function(col_id) {
                                          df[row_id, col_id]
                                        })

                       names(values) <- names(df)[apply_cols]
                     })

  if(!is.null(name_col)) {
    names(out_list) <- df[[name_col]]
  }

  out_list
}
