#' Quality-control metrics calculation
#'
#' Calculate QC metrics, includin Z' score, CV of DMSO negative control, per assay plate.
#'
#' @param d data.frame contains essential columns with the raw data.
#' @param resp response type, specify either 'nval_median' or 'nval_mean' for QC calculation
#' @param assay_info assay_info list, contains names of primary and cytotox assay, names must match
#' what are provided in the raw data, under the column 'assay'.
#'
#' @examples
#' ## calculate QC measures from demo data
#' assay_info <- list(prim_assay = "Primary",toxi_assay = "Cytotox")
#' demo_mc_norm <- normalize_per_plate(demo_mc, nctrl = "DMSO")
#' qc <- qc_per_plate(demo_mc_norm, assay_info)
#'
#' @return three dataframe each representing negative control stats, positive control stats and QC metrics (CV_DMSO, Z' score, SSMD) for each assay plate
#' @export
# QC measures are calculated per 96-well palte. 'apid' plus 'assay' column can serve as the ID to distinguish each plate.

qc_per_plate <- function(d, assay_info, resp = "nval_median") {

  # claim variables for passing R CMD Check
  assay <- apid <- wllt <- conc <- CV_DMSO <- SSMD <- Z_prime <- logc <- mean_DMSO <- mean_positive <- NULL
  sd_positive <- sd_DMSO <- unique_id <- NULL



  # choose which normalized value to use (median based or mean based)
  if (resp == "nval_mean") {
    d$resp <- d$nval_mean
  } else if (resp == "nval_median") {
    d$resp <- d$nval_median
  } else {
    stop("specify either nval_median or nval_mean for QC calculation")
  }

  # create unique ID for each plate (apid + assay)
  # (can also distinguish between primary and cytotox assays)
  # group by the unique ID.
  d <- d %>% #dplyr::mutate(uid = paste(apid, assay)) %>%
      # dplyr::group_by(pid, assay, repi)
      dplyr::group_by(apid, assay)

  # calculate vehicle control (e.g., DMSO etc.) qc stats
  # use welltype = n to identify vehicle control
  n_ctrl_sum <- d %>%
    dplyr::filter(wllt == "n") %>%
    dplyr::summarize(count_DMSO = n(),                     #count of DMSO control wells on each plate
                     count_DMSO_NA = sum(is.na(resp)),     #count of DMSO control wells with missing values
                     mean_DMSO = mean(resp, na.rm = TRUE),
                     sd_DMSO = stats::sd(resp, na.rm = TRUE),
                     CV_DMSO = 100 * stats::sd(resp, na.rm = TRUE) / mean(resp, na.rm = TRUE)
                     # median_DMSO= median(resp, na.rm=TRUE),
                     # mad_DMSO= mad(resp, constant = 1, na.rm=TRUE),
                     # bmad_DMSO= mad(resp, constant = 1.4826, na.rm=TRUE),
                     # three_bmad = 3*bmad_DMSO
                     )


  # calculate positive control stats
  pctrl_prim <-  d %>%
    dplyr::filter(assay == assay_info$prim_assay,
                  wllt == "pr") %>%
    dplyr::filter(conc == max(conc)) %>%
    dplyr::summarize(sd_positive = stats::sd(resp, na.rm = TRUE),
                     mean_positive= mean(resp, na.rm = TRUE))

  pctrl_toxi <-  d %>%
    dplyr::filter(assay == assay_info$toxi_assay,
                  wllt == "pc") %>%
    dplyr::filter(conc == max(conc)) %>%
    dplyr::summarize(sd_positive = stats::sd(resp, na.rm = TRUE),
                     mean_positive = mean(resp, na.rm = TRUE))

  p_ctrl_sum <- bind_rows(pctrl_prim, pctrl_toxi)
  remove(pctrl_toxi, pctrl_prim)


  # calculate Z', SSMD

  qc <- dplyr::left_join(p_ctrl_sum, n_ctrl_sum,
                        by=c("apid"="apid", "assay"="assay")) %>%
       #replace NaN with 0 in sd_positive
       dplyr::mutate(sd_positive = (ifelse(is.na(sd_positive), 0, sd_positive))) %>%
       dplyr::mutate(Z_prime = 1 - 3*(sd_positive + sd_DMSO) / (abs(mean_positive - mean_DMSO)),
                       SSMD = (abs(mean_positive - mean_DMSO) / sqrt(sd_positive^2 + sd_DMSO^2))) %>%
       dplyr::mutate(unique_id = paste(apid, assay, sep = "_")) %>%
       dplyr::select(unique_id, apid, assay, CV_DMSO, Z_prime, SSMD)



  list(neg_ctrl_sum = n_ctrl_sum,
       pos_ctrl_sum = p_ctrl_sum,
       qc = qc)
}


