#' normalize per plate
#'
#' normalize raw readings as percent of median vehicle control wells
#'
#' @param dt data.frame contains essential columns with the raw data.
#' @param nctrl the name (spid) of the vehicle/solvent control used for calculation
#'
#' @return data.frame with normalized value columns. 'nval_mean' column is the normalized value
#' calculated using the mean of vehicle control wells, 'nval_median' column is the normalized value
#' calculated using the median of vehicle control wells.
#' @import dplyr
#'
#' @examples
#' ## normalize demo data
#' demo_mc_norm <- normalize_per_plate(demo_mc, nctrl = "DMSO")
#'
#' @export

normalize_per_plate <- function(dt, nctrl = "DMSO") {

  # claim variables for passing R CMD Check
  assay <- spid <- apid <- rval <- NULL


  output <- data.frame()
  #iterate through each assay, calculate separately
  for (as in unique(dt$assay)) {
    df_t <- normalize_single_assay(dplyr::filter(dt, assay == as), nctrl)
    output <- dplyr::bind_rows(output, df_t)
  }

  return(output)

}

#' normalize per plate (from single assay data)
#'
#' normalize raw readings as percent of median/mean vehicle control wells, per assay, per plate.
#' This function is called by normalize_per_plate, should not be called directly by user.
#' #'
#' @param dt data.frame contains essential columns with the raw data.
#' @param nctrl the name (spid) of the vehicle/solvent control used for calculation
#'
#' @return data.frame with normalized value columns. 'nval_mean' column is the normalized value
#' calculated using the mean of vehicle control wells, 'nval_median' column is the normalized value
#' calculated using the median of vehicle control wells.
#' @import dplyr
#
normalize_single_assay <- function(dt, nctrl) {

  # claim variables for passing R CMD Check
  assay <- spid <- apid <- rval <- NULL


  #calculate the mean and median value for negative control
  ctrl_avg <- dt %>%
    dplyr::filter(spid == nctrl) %>%
    dplyr::group_by(apid) %>%
    dplyr::summarize(mean_DMSO = mean(rval, na.rm = TRUE),
                     median_DMSO = stats::median(rval, na.rm = TRUE))

  #iterate through the data table to calculate normalized percent activity
  #normalization was done using both mean and median separately.

  temp <- data.frame()
  for (id in ctrl_avg$apid){
    med <- dplyr::filter(ctrl_avg, apid == id)$median_DMSO
    avg <- dplyr::filter(ctrl_avg, apid == id)$mean_DMSO
    t <- dplyr::filter(dt, apid == id) %>%
      mutate(nval_mean = 100 * rval/avg,
             nval_median = 100 * rval/med)
    temp <- bind_rows(temp, t)
  }
  return(temp)
}
