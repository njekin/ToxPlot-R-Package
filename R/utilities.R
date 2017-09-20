
#' save plots in pdf
#'
#' save ggplot2 plots generated in a list to a pdf file
#'
#' @param plot_list  the r list object contains all ggplot2 objects
#' @param filename  the output file name, including the file directory
#'
#' @examples
#' ## start from raw data
#' # define assay
#' assay_info <- list(prim_assay = "Primary",toxi_assay = "Cytotox")
#' # data normalization
#' demo_mc_norm <- normalize_per_plate(demo_mc, nctrl = "DMSO")
#' # filter out two test chemicals
#' demo_mc_norm <- dplyr::filter(demo_mc_norm, spid %in% c("TP0001502B05", "TP0001502B01"))
#' # fit curve with default 20% threshold
#' demo_md <- fit_curve_tcpl(demo_mc_norm, assay_info)
#' # calculate TAA and Med_diff only
#' demo_rank <- rank_tcpl(demo_md, med_taa = NULL, med_med_diff = NULL)
#' #produce plots with notations
#' demo_plots <- plot_tcpl_minimal(demo_md, demo_rank, notation = TRUE)
#'
#' ## save all the plots as pdf
#' # save_plot_pdf(demo_plots, ".\output plots\all_plots.pdf")
#'
#' ## save the 1st plot as pdf
#' # save_plot_pdf(demo_plots[1], ".\output plots\plot1.pdf")
#'
#'
#'
#' @export
#'

save_plot_pdf <- function(plot_list, filename) {
  # claim variables for passing R CMD Check
  assay <- spid <- apid <- rval <- NULL

  grDevices::pdf(filename)
  cat("Preparing file...\n")
  invisible(lapply(plot_list, print))
  grDevices::dev.off()
  cat("Finished!")
}



#' round digits of numbers
#'
#' round numbers in a datafram to specified digits
#'
#' @param df the dataframe input
#' @param digits  the specified number of digits
#' @return a dataframe
#'
#'
round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  df[,nums] <- round(df[,nums], digits = digits)
  df
}


# -----------------------------------------------------
# function to scale dataset from 0 to 100
# -----------------------------------------------------

rescale_0_100 <- function(x) {
  y <- scale(x, center = min(x, na.rm = T), scale = diff(range(x, na.rm = T)))
  as.vector(y) * 100
}



## function to generate index for concentrations

## both are targeting unique(x), 'match' provide the position of each x value on the list of unique(x),
## while 'rank' provide the ordering index that is exported for each component of unique(x)

generate_index <- function(x) {
  uls <- unique(x)
  index <- rank(uls)[match(x, uls)]
  return(as.integer(index))
}
