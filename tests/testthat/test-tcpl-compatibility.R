test_that("tcpl results preserve scalar, vector and absent predictions", {
  fit <- list(hill = NA_integer_, logc = 1:18, hill_modl = NULL,
              gnls_modl = NULL, hill_tp = NA_real_)
  result <- toxplot:::tcpl_fit_row(fit)
  expect_equal(nrow(result), 1L)
  expect_identical(result$logc[[1]], 1:18)
  expect_null(result$hill_modl[[1]])
  expect_null(result$gnls_modl[[1]])
  expect_identical(result$hill, NA_integer_)
})

test_that("the complete archived failing example fits, ranks and renders", {
  info <- list(prim_assay = "Primary", toxi_assay = "Cytotox")
  models <- fit_curve_tcpl(demo_mc_norm, info)
  expect_length(models, length(unique(demo_mc_norm$spid)))
  for (m in models) {
    for (assay in c("prim", "toxi")) {
      fit <- m[[paste0("model_", assay)]]
      dat <- m[[paste0("data_", assay)]]
      expect_equal(nrow(fit), 1L)
      expect_identical(fit$logc[[1]], dat$logc)
      expect_equal(fit$apid, dat$apid[1])
    }
  }
  ranked <- rank_tcpl(models)
  expect_equal(nrow(ranked), length(models))
  expect_equal(nrow(summary_tcpl(models)), length(models))
  grDevices::pdf(file = NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (plotter in list(plot_tcpl, plot_tcpl_minimal)) {
    plots <- plotter(models, ranked, notation = TRUE)
    expect_length(plots, length(models))
    expect_warning(lapply(plots, ggplot2::ggplotGrob), NA)
  }
})

test_that("inactive and insufficient-concentration series remain one row", {
  info <- list(prim_assay = "Primary", toxi_assay = "Cytotox")
  d <- subset(demo_mc_norm, spid == "TP0001502B03")
  d$nval_median <- 100
  for (dat in list(d, subset(d, conc == min(conc)))) {
    result <- fit_curve_tcpl(dat, info)
    expect_equal(nrow(result[[1]]$model_prim), 1L)
    expect_true(is.na(result[[1]]$model_prim$hill))
    expect_null(result[[1]]$model_prim$hill_modl[[1]])
    expect_equal(nrow(rank_tcpl(result)), 1L)
  }
})

test_that("PDF export closes its device even if printing fails", {
  before <- grDevices::dev.cur()
  output <- tempfile(fileext = ".pdf")
  on.exit(unlink(output), add = TRUE)
  save_plot_pdf(list(ggplot2::ggplot()), output)
  expect_identical(grDevices::dev.cur(), before)
  expect_gt(file.info(output)$size, 0)
  bad <- ggplot2::ggplot(data.frame(x = 1),
                         ggplot2::aes(x = x, y = not_a_column)) +
    ggplot2::geom_point()
  expect_error(save_plot_pdf(list(bad), output), "not_a_column")
  expect_identical(grDevices::dev.cur(), before)
})
