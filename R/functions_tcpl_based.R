#' hill mode in ToxCast tcpl package
#'
#' @param p a vector containing the Hill model parameters: top, log AC50, hill coefficient
#' @param x a vector of log concentrations
#' @return calculated y value based on the x and model parameters
#'
hill_model <- function(p, x) {
  y <- p[1] / (1 + 10 ^ ((p[2] - x) * p[3]))
  return(y)
}

# ##========inversed hill model (using 100- the Y value)=======##
# hill_model_inverse <- function(p, x) {
#   ###   p:     a numeric vector of length 4 containg the starting values for
#   ###          the hill model, in order: top, log AC50, hill
#   ###          coefficient
#   ###   x: a numeric vector containing the log concentration values
#
#   y <- 100 - (p[1] / (1 + 10 ^ ((p[2] - x) * p[3])))
#   return(y)
# }

#' calculate absolute EC_anything based on tcpl hill model
#'
#' @param p a vector containing the Hill model parameters: top, log AC50, hill coefficient
#' @param y the y value
#' @return calculated x value

log_abs_ec <- function(p, y) {
  ## y value is the absolute response value
  ##p[1] is top parameter
  ##p[2] is ga, logac50
  ##p[3] is gw, hillslope
  suppressWarnings(x <- p[2] - log10((p[1] - y) / y) / p[3])
  if (is.nan(x)) {x <- NA}
  return(x)
}


#' function to calculate Area Under the Curve (AUC) of the hill model
#'
#' @param p a vector containing the Hill model parameters: top, log AC50, hill coefficient
#' @param lower lower boundary of x for integration
#' @param upper upper boundary of x for integration
#' @return calculated area under the curve (AUC) value
#'
auc_hill_tcpl <- function(p, lower, upper) {
  #define the hill model function
  hill <- function(x) {
    p[1] / (1 + 10 ^ ((p[2] - x) * p[3]))
  }
  a <- stats::integrate(hill, lower = lower, upper = upper)
  return(a[[1]])
}




#' fit dose-resopnse curve using tcpl hill model
#'
#' Curve fitting using the tcplFit function in `tcpl` package.
#' Chemicals are modelled based on spid.
#' If you want to model the same chemical (e.g. positive controls),
#' then assign different spid to this chemical so the function can separate them out.
#' Absolute IC20 and absolute IC50 are calculated as well.
#'
#' @param df input data contain normalized assay readings
#' @param assay_info predefined names for primary and cytotoxicity assays,
#' use NULL if either one of the assay does not need to be modeled.
#' @param prim_cutoff significance cutoff for primary assay (eg. 3sigma or 3bMAD)
#' @param toxi_cutoff significance cutoff for cytotoxicity assay (eg. 3sigma or 3bMAD)
#'
#' @return A list object containing modeling results, the corresponding data for each chemical.
#'
#' @examples
#' ## fit curve with default significant threshold 20
#'
#' demo_md <- fit_curve_tcpl(mc_norm, assay_info =
#' list(prim_assay = "Primary", toxi_assay = "Cytotox"))
#'
#' ## start from raw data
#' # define assay
#' assay_info <- list(prim_assay = "Primary",toxi_assay = "Cytotox")
#' # data normalization
#' demo_mc_norm <- normalize_per_plate(demo_mc, nctrl = "DMSO")
#' # filter out two test chemicals
#' demo_mc_norm <- dplyr::filter(demo_mc_norm, spid %in% c("TP0001502B05", "TP0001502B01"))
#' # fit curve with default 20% threshold
#' demo_md <- fit_curve_tcpl(demo_mc_norm, assay_info)
#'
#' ## fit curve with specified significance threshold
#' demo_md <- fit_curve_tcpl(demo_mc_norm, assay_info, prim_cutoff = 25, toxi_cutoff = 25)
#'
#'
#' @import dplyr
#' @export

fit_curve_tcpl <- function(df, assay_info, prim_cutoff = 20, toxi_cutoff = 20) {

  # claim variables for passing R CMD Check
  conc <- spid <- assay <- apid <- NULL

  st_time <- Sys.time()
  r_list <- list()
  c_list <- list()
  bmad_prim <- prim_cutoff/3
  bmad_toxi <- toxi_cutoff/3
  # get log concentration
  df <- df %>% dplyr::mutate(logc = log10(conc))
  spid_list <- unique(df$spid)
  if (is.null(assay_info$prim) & is.null(assay_info$toxi_assay)) {
    stop("assay_info cannot be NULL for both primary and cytotoxicity assay")
  }
  n <- 1
  prim_md <- toxi_md <- model_list <- list()
  cat("Processing", length(unique(df$spid)), "samples(spid)....\n")
  # process by each spid
  for (id in unique(df$spid)) {
    #d <- df %>% dplyr::filter(spid == i)
    cat(id, "||")
    # model raiu data
    # check if toxi assay data is available, if not, skip modeling
    if (is.null(assay_info$prim_assay)) {
      prim_md <- NA
      prim_dt <- NA
    } else {

      prim_dt <- df %>% dplyr::filter(spid == id, assay == assay_info$prim_assay)
      m <- tcpl::tcplFit(logc = prim_dt$logc, resp = 100 - prim_dt$nval_median, bmad_prim)
      absIC50 <- log_abs_ec(c(m$hill_tp, m$hill_ga, m$hill_gw), 50)
      absIC20 <- log_abs_ec(c(m$hill_tp, m$hill_ga, m$hill_gw), 20)
      m[["absIC20"]] <- absIC20
      m[["absIC50"]] <- absIC50
      m[["apid"]] <- prim_dt[[1,1]]
      m[["assay"]] <- assay_info$prim_assay
      m[["spid"]] <- id
      #prim_md <- dplyr::bind_rows(prim_md, m)
      # print(m)
      prim_md <- data.frame(m) %>% dplyr::select(apid, assay, spid, everything())
    }

    # model cytotox data
    # check if toxi assay data is available, if not, skip modeling
    if (is.null(assay_info$toxi_assay)) {
      toxi_md <- NA
      toxi_dt <- NA
    } else {
      toxi_dt <- df %>% dplyr::filter(spid == id, assay==assay_info$toxi_assay)
      m <- tcpl::tcplFit(logc = toxi_dt$logc, resp = 100 - toxi_dt$nval_median, bmad_toxi)
      absIC50 <- log_abs_ec(c(m$hill_tp, m$hill_ga, m$hill_gw), 50)
      absIC20 <- log_abs_ec(c(m$hill_tp, m$hill_ga, m$hill_gw), 20)
      m[["absIC20"]] <- absIC20
      m[["absIC50"]] <- absIC50
      m[["apid"]] <- toxi_dt[[1,1]]
      m[["assay"]] <- assay_info$toxi_assay
      m[["spid"]] <- id
      #toxi_md <- dplyr::bind_rows(toxi_md, m)

      toxi_md <- data.frame(m) %>% dplyr::select(apid, assay, spid, everything())

    }

    #build final model list
    model_list[[n]] <- list(spid = id,
                            model_prim = prim_md,
                            model_toxi = toxi_md,
                            data_prim = prim_dt,
                            data_toxi = toxi_dt,
                            cutoff_prim = prim_cutoff,
                            cutoff_toxi = toxi_cutoff,
                            assay_info = assay_info)

    n <- n + 1

  }

  time <- difftime(Sys.time(), st_time) %>% round(1)
  cat("\nCurve Fitting Completed!\nCalculation time:", paste(unclass(time), units(time)), "\n\n")


  return(model_list)
}



#' function to calculate ranking score, TAA, med_diff, EC values based on tcpl hill model
#'
#' calculate ranking score, TAA, med_diff, absolute EC values, AC50, based on the hill model in tcpl package
#'
#' @param tcpl_models the list object returned by 'fit_curve_tcpl' function
#' @param spid_chnm_table a reference table with 'spid' and the corresponding chemical name 'chnm' column,
#' and the CAS number 'casn' column.
#' @param med_taa the median TAA value from reference chemical, if not supplied, then ranking score won't be calculated.
#' @param med_med_diff the median Median-Difference from reference chemical, if not supplied, then ranking score won't be calculated.
#' @return a dataframe containing ranking metrics for each chemical (spid)
#'
#' @examples
#' ## start with normalized data
#' demo_md <- fit_curve_tcpl(mc_norm, assay_info =
#' list(prim_assay = "Primary", toxi_assay = "Cytotox"))
#' demo_rank <- rank_tcpl(demo_md)
#'
#'
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
#'
#' ## calculate ranking score with specified median TAA and meidan Med_Difference
#' demo_rank <- rank_tcpl(demo_md, med_taa = 150, med_med_diff = 92)
#'
#' @export
#'
rank_tcpl <- function(tcpl_models, spid_chnm_table = NULL, med_taa = NULL, med_med_diff = NULL) {

  # claim variables for passing R CMD Check
  taa_norm <- med_diff_norm <- ranking_score <- spid <- logc <- nval_median <- NULL

  df <- data.frame()
  ##iterate through each spid
  for (i in seq_along(tcpl_models)) {
    #get chemical sample id 'spid'
    .spid <- tcpl_models[[i]]$spid
    # print(i)
    # print(.spid)
    #get chnm
    if (!is.null(spid_chnm_table)) {
      chnm <- dplyr::filter(spid_chnm_table, spid == .spid)$chnm
      casn <- dplyr::filter(spid_chnm_table, spid == .spid)$casn
    } else {
      chnm <- NA
      casn <- NA
    }

    assay_info <- tcpl_models[[i]]$assay_info
    if (is.null(assay_info$toxi_assay)) {
      stop("No toxi_assay, ranking score is not calculable!")
    }
    #extract modelling results
    m_toxi <- tcpl_models[[i]]$model_toxi
    m_prim <- tcpl_models[[i]]$model_prim

    #cutoff
    cutoff_prim <- tcpl_models[[i]]$cutoff_prim
    cutoff_toxi <- tcpl_models[[i]]$cutoff_toxi

    #
    #calculate med_diff at max(logc)
    max_prim <- tcpl_models[[i]]$data_prim %>%
      dplyr::filter(logc == max(logc)) %>%
      dplyr::summarize(median = stats::median(nval_median))
    max_toxi <- tcpl_models[[i]]$data_toxi %>%
      dplyr::filter(logc == max(logc)) %>%
      dplyr::summarize(median = stats::median(nval_median))
    med_diff <- max_toxi$median - max_prim$median

    ##get auc, logEC_3bmadprim, for prim
    if (!is.na(m_prim$hill) & m_prim$hill == 1) {
      para_prim <- c(m_prim$hill_tp, m_prim$hill_ga, m_prim$hill_gw)
      #get auc
      lr <- log_abs_ec(para_prim, cutoff_prim)
      aa_prim <-
        auc_hill_tcpl(para_prim, lower = lr, upper = m_prim$logc_max) - cutoff_prim * (m_prim$logc_max - lr)
      #get AC50
      AC50_prim <- m_prim$hill_ga
      #get absEC
      if (hill_model(para_prim,m_prim$logc_max) > 50) {
        absEC50_prim <- log_abs_ec(para_prim, 50)
        absEC80_prim <- log_abs_ec(para_prim, 20)
      } else if (hill_model(para_prim,m_prim$logc_max) > 20) {
        absEC50_prim <- NA
        absEC80_prim <- log_abs_ec(para_prim, 20)
      } else {
        absEC80_prim <- NA
        absEC50_prim <- NA
      }
    } else {
      aa_prim <- NA
      lr <- NA
      absEC80_prim <- NA
      absEC50_prim <- NA
      AC50_prim <- NA
    }

    ##get auc, logEC_3bmadprim, for toxitox
    if (!is.na(m_toxi$hill) & m_toxi$hill == 1) {
      para_toxi <- c(m_toxi$hill_tp, m_toxi$hill_ga, m_toxi$hill_gw)
      #get AC50
      AC50_toxi <- m_toxi$hill_ga
      #get aa_toxi
      if (hill_model(para_toxi,m_prim$logc_max) > cutoff_prim) {
        lc <- log_abs_ec(para_toxi, cutoff_prim)
        aa_toxi <-
          auc_hill_tcpl(para_toxi,lower = lc, upper = m_prim$logc_max) - cutoff_prim * (m_prim$logc_max - lc)
      } else {
        aa_toxi <- NA
        lc <- -3
      }

      #calculate absEC
      if (hill_model(para_toxi,m_prim$logc_max) > 50) {
        absEC50_toxi <- log_abs_ec(para_toxi, 50)
        absEC80_toxi <- log_abs_ec(para_toxi, 20)
      } else if (hill_model(para_toxi,m_prim$logc_max) > 20) {
        absEC50_toxi <- NA
        absEC80_toxi <- log_abs_ec(para_toxi, 20)
      } else {
        absEC80_toxi <- NA
        absEC50_toxi <- NA
      }
      #get cytotox limit (absEC_cutoff_toxi)
      absEC_ct_toxi <- log_abs_ec(para_toxi, cutoff_toxi)
    } else {
      aa_toxi <- NA
      lc <- -3
      absEC80_toxi <- NA
      absEC50_toxi <- NA
      AC50_toxi <- NA
      absEC_ct_toxi <- NA
    }

    ##calculate taa
    if (is.na(aa_prim)) {
      taa <- NA
    } else if (is.na(aa_toxi)) {
      taa <- aa_prim
    } else {
      taa <- aa_prim - aa_toxi
    }

    ##calculate selectivity based logEC value at 3bmad of prim assay.
    if (is.na(lr)) {
      selectivity_3bmad <- NA
    } else {
      selectivity_3bmad <- lc - lr
    }

    ##calculate selectivity based on logAC50 (original method)
    if (is.na(m_prim$hill_ga)) {
      selectivity_AC50 <- NA
    } else if (is.na(m_toxi$hill_ga)) {
      selectivity_AC50 <- -3 - m_prim$hill_ga
    } else {
      selectivity_AC50 <- m_toxi$hill_ga - m_prim$hill_ga
    }


    ##gather taa and selectivity into one table
    t <- data.frame(
      index = i,
      spid = .spid,
      chnm = chnm,
      casn = casn,
      # aa_toxi = aa_toxi,
      # aa_prim = aa_prim,
      taa = taa,
      # sel_3bMAD = selectivity_3bmad,
      # sel_AC50 = selectivity_AC50,
      med_diff = med_diff,
      AC50_toxi = AC50_toxi,
      AC50_prim = AC50_prim,
      absEC80_toxi = absEC80_toxi,
      absEC50_toxi = absEC50_toxi,
      absEC80_prim = absEC80_prim,
      absEC50_prim = absEC50_prim,
      cyto_lim = absEC_ct_toxi
    )

    df <- base::rbind(df, t)

  }

  #calculate ranking_score
  #this done by adding 0-100 rescaled TAA value an med_diff value
  # df <- df %>%
  #   dplyr::mutate(taa_rescale = rescale_0_100(taa),
  #                 med_diff_rescale = rescale_0_100(med_diff),
  #                 ranking_score = taa_rescale + med_diff_rescale ) %>%
  #   dplyr::arrange(desc(ranking_score))
  if (is.null(med_taa)&is.null(med_med_diff)) {
    df <- df %>%
      dplyr::mutate(ranking_score = NA)
  } else {
    df <- df %>%
      dplyr::mutate(taa_norm = taa/med_taa*100,
                    med_diff_norm = med_diff/med_med_diff*100,
                    ranking_score = taa_norm + med_diff_norm ) %>%
      dplyr::arrange(desc(ranking_score))

  }



  return(df)
}




#' function to summarize curve fitting results
#'
#' @param tcpl_models the list object returned by 'fit_curve_tcpl' function
#' @param spid_chnm_table a reference table with 'spid' and the corresponding chemical name 'chnm' column,
#' and the CAS number 'casn' column.
#' @return a data.frame contains summarized metrics for each chemical (spid)
#'
#' @examples
#' ## supply models as the essential argument. spid_chnm_table is optional.
#' demo_md <- fit_curve_tcpl(mc_norm, assay_info =
#' list(prim_assay = "Primary", toxi_assay = "Cytotox"))
#' demo_sum <- summary_tcpl(demo_md)
#'
#'
#' ## start from raw data
#' # define assay
#' assay_info <- list(prim_assay = "Primary",toxi_assay = "Cytotox")
#' # data normalization
#' demo_mc_norm <- normalize_per_plate(demo_mc, nctrl = "DMSO")
#' # filter out two test chemicals
#' demo_mc_norm <- dplyr::filter(demo_mc_norm, spid %in% c("TP0001502B05", "TP0001502B01"))
#' # fit curve with default 20% threshold
#' demo_md <- fit_curve_tcpl(demo_mc_norm, assay_info)
#' # obtain summary table
#' demo_sum <- summary_tcpl(demo_md)
#'
#' @export
#'
summary_tcpl <- function(tcpl_models, spid_chnm_table = NULL) {

  # claim variables for passing R CMD Check
  spid <- NULL


  df <- data.frame()
  ##iterate through each spid
  for (i in seq_along(tcpl_models)) {
    #get chemical sample id 'spid'
    .spid <- tcpl_models[[i]]$spid
    # print(i)
    # print(.spid)
    #get chnm
    if (!is.null(spid_chnm_table)) {
      chnm <- dplyr::filter(spid_chnm_table, spid == .spid)$chnm
      casn <- dplyr::filter(spid_chnm_table, spid == .spid)$chnm
    } else {
      chnm <- NA
    }

    #extract modelling results
    m_toxi <- tcpl_models[[i]]$model_toxi
    m_prim <- tcpl_models[[i]]$model_prim

    #cutoff
    cutoff_prim <- tcpl_models[[i]]$cutoff_prim
    cutoff_toxi <- tcpl_models[[i]]$cutoff_toxi

    #get assay_info
    assay_info <- tcpl_models[[i]]$assay_info

    ##get metrics for prim
    if (is.null(assay_info$prim_assay)) {
      absEC80_prim <- NA
      absEC50_prim <- NA
      AC50_prim <- NA
    } else {
      if (!is.na(m_prim$hill) & m_prim$hill == 1) {
        para_prim <- c(m_prim$hill_tp, m_prim$hill_ga, m_prim$hill_gw)
        #get AC50
        AC50_prim <- m_prim$hill_ga
        #get absEC
        if (hill_model(para_prim, m_prim$logc_max) > 50) {
          absEC50_prim <- log_abs_ec(para_prim, 50)
          absEC80_prim <- log_abs_ec(para_prim, 20)
        } else if (hill_model(para_prim,m_prim$logc_max) > 20) {
          absEC50_prim <- NA
          absEC80_prim <- log_abs_ec(para_prim, 20)
        } else {
          absEC80_prim <- NA
          absEC50_prim <- NA
        }
      } else {
        absEC80_prim <- NA
        absEC50_prim <- NA
        AC50_prim <- NA
      }
    }


    ##get metrics for toxi

    if (is.null(assay_info$toxi_assay)) {
      absEC80_toxi <- NA
      absEC50_toxi <- NA
      AC50_toxi <- NA
    } else {
      if (!is.na(m_toxi$hill) & m_toxi$hill == 1) {
        para_toxi <- c(m_toxi$hill_tp, m_toxi$hill_ga, m_toxi$hill_gw)
        #get AC50
        AC50_toxi <- m_toxi$hill_ga
        #calculate absEC
        if (hill_model(para_toxi,m_toxi$logc_max) > 50) {
          absEC50_toxi <- log_abs_ec(para_toxi, 50)
          absEC80_toxi <- log_abs_ec(para_toxi, 20)
        } else if (hill_model(para_toxi,m_toxi$logc_max) > 20) {
          absEC50_toxi <- NA
          absEC80_toxi <- log_abs_ec(para_toxi, 20)
        } else {
          absEC80_toxi <- NA
          absEC50_toxi <- NA
        }
      } else {
        absEC80_toxi <- NA
        absEC50_toxi <- NA
        AC50_toxi <- NA
      }

    }




    ##gather taa and selectivity into one table
    t <- data.frame(
      index = i,
      spid = .spid,
      chnm = chnm,
      AC50_toxi = AC50_toxi,
      AC50_prim = AC50_prim,
      absEC80_toxi = absEC80_toxi,
      absEC50_toxi = absEC50_toxi,
      absEC80_prim = absEC80_prim,
      absEC50_prim = absEC50_prim
    )

    df <- base::rbind(df, t)

  }


  return(df)
}




#' Plot dose-resonse curves based on the tcpl hill model
#'
#' Produce the plot for the dose-response curves and data points for both primary and toxicity assay.
#' The direction of the data and dose-resonse curves are presented as the original data, rather than
#' the uptrend direction required by the 'tcpl' function. Plots are sorted by the ranking_score.
#'
#' @param tcpl_models the list object created by 'fit_curve_tcpl' function
#' @param rank_table the data.frame output from 'rank_tcpl' function
#' @param spid_chnm_table the spid, chnm, casn info table
#' @param notation value can be TRUE or FALSE, determine whehter to show potency metrics on the plot
#' @param cunit  the unit of concentration, on default is "M" (molar).
#'
#' @return list of ggplot2 objects, each corresponding to one spid.
#' @import ggthemes ggplot2
#'
#' @examples
#' ## produce plots without notations
#' demo_md <- fit_curve_tcpl(mc_norm, assay_info =
#' list(prim_assay = "Primary", toxi_assay = "Cytotox"))
#' plots <- plot_tcpl(demo_md)
#'
#'
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
#' demo_plots <- plot_tcpl(demo_md, demo_rank, notation = TRUE)
#'
#' ##produce plots with notations, with changed concentration unit displayed on the plot
#' demo_plots <- plot_tcpl(demo_md, demo_rank, notation = TRUE, cunit = "uM")
#'
#' @export
#'
#'
plot_tcpl <-
  function(tcpl_models, rank_table=NULL, spid_chnm_table = NULL, notation = FALSE, cunit = "M") {


    # claim variables for passing R CMD Check
    spid <- nval_median <- logc <- resp <- pred <- assay <- NULL


    #initiate empty output plot list
    plot_list <- list()
    #reorder if rank is available
    if (!is.null(rank_table)) {
      tcpl_models <- tcpl_models[rank_table$index]
    }

    #loop through tcpl_models's unique spid.
    #note tcpl_models is the level 4 output from tcpl package, including all the modelling results
    for (i in seq_along(tcpl_models)) {
      #get chemical sample id and chemical name
      .spid <- tcpl_models[[i]]$spid
      #print(.spid)
      if (!is.null(spid_chnm_table)) {
        chnm <- dplyr::filter(spid_chnm_table, spid == .spid)$chnm
        casn <- dplyr::filter(spid_chnm_table, spid == .spid)$casn
      } else {
        chnm <- NA
        casn <- NA
      }


      #cutoff
      cutoff_prim <- tcpl_models[[i]]$cutoff_prim
      cutoff_toxi <- tcpl_models[[i]]$cutoff_toxi
      #assay_info
      assay_info <- tcpl_models[[i]]$assay_info

      #get normalized resonse value (cytotox and raiu together)
      #and use 100 minus the response value (inverse the plot)
      d <- dplyr::bind_rows(tcpl_models[[i]]$data_prim, tcpl_models[[i]]$data_toxi) %>%
        dplyr::mutate(resp = nval_median) #%>%
      #mutate(assay = ifelse(aeid == 1, assay_info$toxi_assay, assay_info$prim_assay))


      #determine left and right x boundary of the plot
      if ( round(max(d$logc)) < max(d$logc) ) {
        rb <- round(max(d$logc)) + 1
      } else { rb <- round(max(d$logc))}

      lb <- round(min(d$logc))-1

      #initiate basic plot with data points.
      g <- ggplot(d, aes(x = logc, y = resp))

      if (!is.null(spid_chnm_table)){
        g <- g + labs(
          title = paste(i, ". SPID: ", .spid, "\nNAME: ", chnm, "\nCAS NO: ", casn, sep = ""),
          x = paste("Concentration (log", cunit, ")", sep = ""),
          y = "% Control Activity"
          )
      } else {
        g <- g + labs(
          title = paste(i, ". SPID:" , .spid),
          x = paste("Concentration (log", cunit, ")", sep = ""),
          y = "% Control Activity"
        )
      }

        # labs(
        #   title = paste(i, "\nSPID: " , .spid, "\nNAME: ", chnm, "\nCAS NO: ", casn, sep = ""),
        #   x = "Concentration (logM)",
        #   y = "% Control Activity"
        # )


      #extract modelling results
      m_toxi <- tcpl_models[[i]]$model_toxi
      m_prim <- tcpl_models[[i]]$model_prim


      #create 100 concentrations
      s <- expand.grid(logc = seq(lb, rb, length = 130))

      #test and plot cytotox model
      if (!is.na(m_toxi$hill) & m_toxi$hill == 1) {
        para_cyto <- c(m_toxi$hill_tp, m_toxi$hill_ga, m_toxi$hill_gw)
        p1 <- 100 - hill_model(para_cyto, s)
        p1 <- dplyr::bind_cols(s, data.frame(p1))
        names(p1) <- c("logc", "pred")
        g <- g + geom_line(
          data = p1,
          aes(x = logc, y = pred),
          size = 2,
          alpha = 0.9,
          color = "#e02929"
        )
        #plot the vertical line at 3bmad cutoff
        #g <- g+ geom_vline(xintercept = log_abs_ec(para_cyto, 3*m_toxi$bmad) )
      }


      #test and plot raiu model
      if (!is.na(m_prim$hill) & m_prim$hill == 1) {
        para_raiu <- c(m_prim$hill_tp, m_prim$hill_ga, m_prim$hill_gw)

        p2 <- 100 - hill_model(para_raiu, s)
        p2 <- bind_cols(s, data.frame(p2))
        names(p2) <- c("logc", "pred")
        g <- g + geom_line(
          data = p2,
          aes(x = logc, y = pred),
          size = 2,
          alpha = 0.9,
          color = "#377eb8"
        )
      }


      #draw 3bmad cutoff line for cytotox and raiu respectively
      g <- g +
        geom_hline(
          yintercept = 100 - cutoff_toxi,
          alpha = 0.5,
          size = 0.5,
          linetype = "dashed",
          color = "#e02929"
        )

      g <- g +
        geom_hline(
          yintercept = 100 - cutoff_prim,
          alpha = 0.5,
          size = 0.5,
          linetype = "dashed",
          color = "#377eb8"
        )

      #plot data points
      #aesthetics fixes
      g <- g +
        geom_point(
          aes(color = assay),
          shape = 21,
          alpha = 0.9,
          size = 3
        ) +
        coord_fixed(
          ylim = c(0, 125),
          #xlim = c(-10, -4),
          ratio = 4 / 120   #used to be 2/70, when x axis was from -9 to -4.
        ) +
        scale_y_continuous(breaks = seq(
          from = 0,
          to = 120,
          by = 20
        )) +
        scale_x_continuous(breaks = seq(
          from = lb,
          to = rb,
          by = 1
        )) +
        theme_few() +
        theme(legend.title = element_blank()) +
        scale_color_manual(values=c("#e02929", "#377eb8"))+
        theme(plot.title=element_text(hjust=0.5))


      ##adding annotations
      if (!is.null(rank_table) & notation == TRUE) {
        #Get ec and ranking info
        ds <- rank_table %>% dplyr::filter(spid==.spid)
        ds <- round_df(ds, digits=2)

        #annotate with text info
        if (is.na(ds$ranking_score)) {
          line1 <- ""} else {
          line1 <- paste("Ranking_Score:", ds$ranking_score)
        }

        line2 <- paste("TAA:", ds$taa)
        line3 <- paste("Med_Diff:", ds$med_diff)
        line4 <- paste(assay_info$prim_assay, "_AC50: ", ds$AC50_prim, sep = "")
        line5 <- paste(assay_info$prim_assay, "_absEC50: ", ds$absEC50_prim, sep = "")
        #line6 <- paste(assay_info$prim_assay, "_absEC80: ", ds$absEC80_prim, sep = "")
        g <- g +
          annotate("text", x= lb + 0.8, y = 30, alpha = 0.8, hjust=0, label=line1) +
          annotate("text", x= lb + 0.8, y = 25, alpha = 0.8, hjust=0, label= line2) +
          annotate("text", x= lb + 0.8, y = 20, alpha = 0.8, hjust=0, label= line3) +
          annotate("text", x= lb + 0.8, y = 15, alpha = 0.8, hjust=0, label= line4) +
          annotate("text", x= lb + 0.8, y = 10, alpha = 0.8, hjust=0, label= line5)
          #annotate("text", x= lb + 0.8, y = 5, alpha = 0.8, hjust=0, label= line6)

      }

      #collect all plots into a list
      plot_list[[i]] <- g


    }

    return(plot_list)

  }



#' Plot dose-resonse curves with minimal text annotation
#' This funciton plots dose-response curve with minimal text annotation,
#' no x and y axis label, 0 borders. Useful when need to present several plots
#' together.
#' @param tcpl_models the list object created by 'fit_curve_tcpl' function
#' @param rank_table the data.frame output from 'rank_tcpl' function
#' @param spid_chnm_table the spid, chnm, casn info table
#' @param notation value can be TRUE or FALSE, determine whehter to show potency metrics on the plot
#' @param cunit  the unit of concentration, on default is "M" (molar).
#'
#' @return list of ggplot2 objects, each corresponding to one spid.
#' @import ggthemes ggplot2
#'
#' @examples
#' ## produce plots without notations
#' demo_md <- fit_curve_tcpl(mc_norm, assay_info =
#' list(prim_assay = "Primary", toxi_assay = "Cytotox"))
#' plots_minimal <- plot_tcpl_minimal(demo_md)
#'
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
#' ##produce plots with notations, with changed concentration unit displayed on the plot
#' demo_plots <- plot_tcpl_minimal(demo_md, demo_rank, notation = TRUE, cunit = "uM")
#'
#'
#'
#'
#' @export
#'
#'
plot_tcpl_minimal <-
  function(tcpl_models, rank_table = NULL, spid_chnm_table = NULL, notation = FALSE, cunit = "M") {


    # claim variables for passing R CMD Check
    spid <- nval_median <- logc <- resp <- pred <- assay <- NULL

    #initiate empty output plot list
    plot_list <- list()

    #reorder if rank is available
    if (!is.null(rank_table)) {
      tcpl_models <- tcpl_models[rank_table$index]
    }

    #loop through tcpl_models's unique spid.
    #note tcpl_models is the level 4 output from tcpl package, including all the modelling results
    for (i in seq_along(tcpl_models)) {
      #get chemical sample id and chemical name
      .spid <- tcpl_models[[i]]$spid
      #print(.spid)
      if (!is.null(spid_chnm_table)) {
        chnm <- dplyr::filter(spid_chnm_table, spid == .spid)$chnm
        casn <- dplyr::filter(spid_chnm_table, spid == .spid)$casn
      } else {
        chnm <- NA
        casn <- NA
      }


      #cutoff
      cutoff_prim <- tcpl_models[[i]]$cutoff_prim
      cutoff_toxi <- tcpl_models[[i]]$cutoff_toxi
      #assay_info
      assay_info <- tcpl_models[[i]]$assay_info

      #get normalized resonse value (cytotox and raiu together)
      #and use 100 minus the response value (inverse the plot)
      d <- dplyr::bind_rows(tcpl_models[[i]]$data_prim, tcpl_models[[i]]$data_toxi) %>%
        dplyr::mutate(resp = nval_median) #%>%
      #mutate(assay = ifelse(aeid == 1, assay_info$toxi_assay, assay_info$prim_assay))


      #determine left and right x boundary of the plot
      if ( round(max(d$logc)) < max(d$logc) ) {
        rb <- round(max(d$logc)) + 1
      } else { rb <- round(max(d$logc))}

      lb <- round(min(d$logc))-1

      #initiate basic plot with data points.
      g <- ggplot(d, aes(x = logc, y = resp))

      if (!is.null(spid_chnm_table)){
        g <- g + labs(
          title = paste(chnm, sep = ""),
          x = paste("Concentration (log", cunit, ")", sep = ""),
          y = "% Control Activity"
        )
      } else {
        g <- g + labs(
          title = paste(i, ". SPID:" , .spid),
          x = paste("Concentration (log", cunit, ")", sep = ""),
          y = "% Control Activity"
        )
      }

      # labs(
      #   title = paste(i, "\nSPID: " , .spid, "\nNAME: ", chnm, "\nCAS NO: ", casn, sep = ""),
      #   x = "Concentration (logM)",
      #   y = "% Control Activity"
      # )


      #extract modelling results
      m_toxi <- tcpl_models[[i]]$model_toxi
      m_prim <- tcpl_models[[i]]$model_prim


      #create 100 concentrations
      s <- expand.grid(logc = seq(lb, rb, length = 130))

      #test and plot cytotox model
      if (!is.na(m_toxi$hill) & m_toxi$hill == 1) {
        para_cyto <- c(m_toxi$hill_tp, m_toxi$hill_ga, m_toxi$hill_gw)
        p1 <- 100 - hill_model(para_cyto, s)
        p1 <- dplyr::bind_cols(s, data.frame(p1))
        names(p1) <- c("logc", "pred")
        g <- g + geom_line(
          data = p1,
          aes(x = logc, y = pred),
          size = 2,
          alpha = 0.9,
          color = "#e02929"
        )
        #plot the vertical line at 3bmad cutoff
        #g <- g+ geom_vline(xintercept = log_abs_ec(para_cyto, 3*m_toxi$bmad) )
      }


      #test and plot raiu model
      if (!is.na(m_prim$hill) & m_prim$hill == 1) {
        para_raiu <- c(m_prim$hill_tp, m_prim$hill_ga, m_prim$hill_gw)

        p2 <- 100 - hill_model(para_raiu, s)
        p2 <- bind_cols(s, data.frame(p2))
        names(p2) <- c("logc", "pred")
        g <- g + geom_line(
          data = p2,
          aes(x = logc, y = pred),
          size = 2,
          alpha = 0.9,
          color = "#377eb8"
        )
      }


      #draw 3bmad cutoff line for cytotox and raiu respectively
      g <- g +
        geom_hline(
          yintercept = 100 - cutoff_toxi,
          alpha = 0.5,
          size = 0.5,
          linetype = "dashed",
          color = "#e02929"
        )

      g <- g +
        geom_hline(
          yintercept = 100 - cutoff_prim,
          alpha = 0.5,
          size = 0.5,
          linetype = "dashed",
          color = "#377eb8"
        )

      #plot data points
      #aesthetics fixes
      g <- g +
        geom_point(
          aes(color = assay),
          shape = 21,
          alpha = 0.9,
          size = 3
        ) +
        coord_fixed(
          ylim = c(0, 125),
          xlim = c(-9, -4),
          ratio = 2 / 70   #used to be 2/70, when x axis was from -9 to -4.
        ) +
        scale_y_continuous(breaks = seq(
          from = 0,
          to = 120,
          by = 20
        )) +
        scale_x_continuous(breaks = seq(
          from = lb,
          to = rb,
          by = 1
        )) +
        theme_few() +
        theme(axis.title.x = element_blank()) +
        theme(axis.title.y = element_blank()) +
        theme(legend.title = element_blank(),
              legend.position = "none") +
        theme(legend.margin=unit(0, "null")) +
        scale_color_manual(values=c("#e02929", "#377eb8"))+
        theme(plot.title=element_text(hjust=0.5, size = 16))


      ##adding annotations
      if (!is.null(rank_table) & notation == TRUE) {
        #Get ec and ranking info
        ds <- rank_table %>% dplyr::filter(spid==.spid)
        ds <- round_df(ds, digits=2)

        #annotate with text info
        line1 <- paste("Ranking_Score:", ds$ranking_score)
        # line2 <- paste("TAA:", ds$taa)
        # line3 <- paste("Med_Diff:", ds$med_diff)
        line4 <- paste("AC50: ", ds$AC50_prim, sep = "")
        line5 <- paste("absEC50: ", ds$absEC50_prim, sep = "")
        # line6 <- paste(assay_info$prim_assay, "_absEC80: ", ds$absEC80_prim, sep = "")
        g <- g +
          annotate("text", x= lb + 0.8, y = 34, alpha = 0.8, hjust=0, label=line1) +
          # annotate("text", x= lb + 0.8, y = 25, alpha = 0.8, hjust=0, label= line2) +
          # annotate("text", x= lb + 0.8, y = 20, alpha = 0.8, hjust=0, label= line3) +
          annotate("text", x= lb + 0.8, y = 22, alpha = 0.8, hjust=0, label= line4) +
          annotate("text", x= lb + 0.8, y = 10, alpha = 0.8, hjust=0, label= line5)
          # annotate("text", x= lb + 0.8, y = 5, alpha = 0.8, hjust=0, label= line6)

      }

      #collect all plots into a list
      plot_list[[i]] <- g


    }

    return(plot_list)

  }


