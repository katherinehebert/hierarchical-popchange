plot_mvgam_fc_custom = function (object, series = 1, newdata, data_test, realisations = FALSE, 
          n_realisations = 15, hide_xlabels = FALSE, xlab, ylab, ylim, 
          n_cores = 1, return_forecasts = FALSE, return_score = FALSE, 
          ...) 
{
  if (!(inherits(object, "mvgam"))) {
    stop("argument \"object\" must be of class \"mvgam\"")
  }
  if (sign(series) != 1) {
    stop("argument \"series\" must be a positive integer", 
         call. = FALSE)
  }
  else {
    if (series%%1 != 0) {
      stop("argument \"series\" must be a positive integer", 
           call. = FALSE)
    }
  }
  if (series > NCOL(object$ytimes)) {
    stop(paste0("object only contains data / predictions for ", 
                NCOL(object$ytimes), " series"), call. = FALSE)
  }
  if (sign(n_realisations) != 1) {
    stop("argument \"n_realisations\" must be a positive integer", 
         call. = FALSE)
  }
  else {
    if (n_realisations%%1 != 0) {
      stop("argument \"n_realisations\" must be a positive integer", 
           call. = FALSE)
    }
  }
  if (return_score) {
    return_forecasts <- TRUE
  }
  if (missing(data_test) & missing("newdata")) {
    if (!is.null(object$test_data)) {
      data_test <- object$test_data
    }
  }
  if (!missing("newdata")) {
    data_test <- newdata
    if (terms(formula(object$call))[[2]] != "y") {
      data_test$y <- data_test[[terms(formula(object$call))[[2]]]]
    }
  }
  data_train <- object$obs_data
  ends <- seq(0, dim(mcmc_chains(object$model_output, "ypred"))[2], 
              length.out = NCOL(object$ytimes) + 1)
  starts <- ends + 1
  starts <- c(1, starts[-c(1, (NCOL(object$ytimes) + 1))])
  ends <- ends[-1]
  if (object$fit_engine == "stan") {
    preds <- mcmc_chains(object$model_output, "ypred")[, 
                                                       seq(series, dim(mcmc_chains(object$model_output, 
                                                                                   "ypred"))[2], by = NCOL(object$ytimes)), drop = FALSE]
  }
  else {
    preds <- mcmc_chains(object$model_output, "ypred")[, 
                                                       starts[series]:ends[series], drop = FALSE]
  }
  s_name <- levels(data_train$series)[series]
  if (!missing(data_test)) {
    if (terms(formula(object$call))[[2]] != "y") {
      if (object$family %in% c("binomial", "beta_binomial")) {
        resp_terms <- as.character(terms(formula(object$call))[[2]])
        resp_terms <- resp_terms[-grepl("cbind", resp_terms)]
        trial_name <- resp_terms[2]
        data_test$y <- data_test[[resp_terms[1]]]
        if (!exists(trial_name, data_test)) {
          stop(paste0("Variable ", trial_name, " not found in newdata"), 
               call. = FALSE)
        }
      }
      else {
        data_test$y <- data_test[[terms(formula(object$call))[[2]]]]
      }
    }
    if (!"y" %in% names(data_test)) {
      data_test$y <- rep(NA, NROW(data_test))
    }
    if (inherits(data_test, "list")) {
      if (!"time" %in% names(data_test)) {
        stop("data_test does not contain a \"time\" column")
      }
      if (!"series" %in% names(data_test)) {
        data_test$series <- factor("series1")
      }
    }
    else {
      if (!"time" %in% colnames(data_test)) {
        stop("data_test does not contain a \"time\" column")
      }
      if (!"series" %in% colnames(data_test)) {
        data_test$series <- factor("series1")
      }
    }
    if (inherits(data_test, "list")) {
      all_obs <- c(data.frame(y = data_train$y, series = data_train$series, 
                              time = data_train$time) %>% dplyr::filter(series == 
                                                                          s_name) %>% dplyr::select(time, y) %>% dplyr::distinct() %>% 
                     dplyr::arrange(time) %>% dplyr::pull(y), data.frame(y = data_test$y, 
                                                                         series = data_test$series, time = data_test$time) %>% 
                     dplyr::filter(series == s_name) %>% dplyr::select(time, 
                                                                       y) %>% dplyr::distinct() %>% dplyr::arrange(time) %>% 
                     dplyr::pull(y))
    }
    else {
      all_obs <- c(data_train %>% dplyr::filter(series == 
                                                  s_name) %>% dplyr::select(time, y) %>% dplyr::distinct() %>% 
                     dplyr::arrange(time) %>% dplyr::pull(y), data_test %>% 
                     dplyr::filter(series == s_name) %>% dplyr::select(time, 
                                                                       y) %>% dplyr::distinct() %>% dplyr::arrange(time) %>% 
                     dplyr::pull(y))
    }
    if (dim(preds)[2] != length(all_obs)) {
      s_name <- levels(object$obs_data$series)[series]
      if (attr(object$model_data, "trend_model") == "None") {
        if (class(object$obs_data)[1] == "list") {
          series_obs <- which(data_test$series == s_name)
          series_test <- lapply(data_test, function(x) {
            if (is.matrix(x)) {
              matrix(x[series_obs, ], ncol = NCOL(x))
            }
            else {
              x[series_obs]
            }
          })
        }
        else {
          series_test = data_test %>% dplyr::filter(series == 
                                                      s_name)
        }
        fc_preds <- predict.mvgam(object, newdata = series_test, 
                                  type = "response", n_cores = n_cores)
      }
      else {
        fc_preds <- forecast.mvgam(object, data_test = data_test, 
                                   n_cores = n_cores)$forecasts[[series]]
      }
      preds <- cbind(preds, fc_preds)
    }
  }
  preds_last <- preds[1, ]
  probs = c(0.05, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.95)
  cred <- sapply(1:NCOL(preds), function(n) quantile(preds[, 
                                                           n], probs = probs, na.rm = TRUE))
  c_light <- c("#DCBCBC")
  c_light_highlight <- c("#C79999")
  c_mid <- c("#B97C7C")
  c_mid_highlight <- c("#A25050")
  c_dark <- c("#8F2727")
  c_dark_highlight <- c("#7C0000")
  if (missing(ylim)) {
    ytrain <- data.frame(series = data_train$series, time = data_train$time, 
                         y = data_train$y) %>% dplyr::filter(series == s_name) %>% 
      dplyr::select(time, y) %>% dplyr::distinct() %>% 
      dplyr::arrange(time) %>% dplyr::pull(y)
    if (tolower(object$family) %in% c("beta", "bernoulli")) {
      ylim <- c(min(cred, min(ytrain, na.rm = TRUE)), 
                max(cred, max(ytrain, na.rm = TRUE)))
      ymin <- max(0, ylim[1])
      ymax <- min(1, ylim[2])
      ylim <- c(ymin, ymax)
    }
    else if (tolower(object$family) %in% c("lognormal", 
                                           "gamma")) {
      ylim <- c(min(cred, min(ytrain, na.rm = TRUE)), 
                max(cred, max(ytrain, na.rm = TRUE)))
      ymin <- max(0, ylim[1])
      ymax <- max(ylim)
      ylim <- c(ymin, ymax)
    }
    else {
      ylim <- c(min(cred, min(ytrain, na.rm = TRUE)), 
                max(cred, max(ytrain, na.rm = TRUE)))
    }
  }
  if (missing(ylab)) {
    ylab <- paste0("Predicitons for ", levels(data_train$series)[series])
  }
  if (missing(xlab)) {
    xlab <- "Time"
  }
  pred_vals <- seq(1:length(preds_last))
  if (hide_xlabels) {
    plot(1, type = "n", bty = "L", xlab = "", xaxt = "n", 
         ylab = ylab, xlim = c(0, length(preds_last)), ylim = ylim, 
         ...)
  }
  else {
    plot(1, type = "n", bty = "L", xlab = xlab, ylab = ylab, 
         xaxt = "n", xlim = c(0, length(preds_last)), ylim = ylim, 
         ...)
    if (!missing(data_test)) {
      axis(side = 1, at = floor(seq(0, max(data_test$time) - 
                                      (min(object$obs_data$time) - 1), length.out = 6)), 
           labels = floor(seq(min(object$obs_data$time), 
                              max(data_test$time), length.out = 6)))
    }
    else {
      axis(side = 1, at = floor(seq(0, max(object$obs_data$time) - 
                                      (min(object$obs_data$time) - 1), length.out = 6)), 
           labels = floor(seq(min(object$obs_data$time), 
                              max(object$obs_data$time), length.out = 6)))
    }
  }
  if (realisations) {
    for (i in 1:n_realisations) {
      lines(x = pred_vals, y = preds[i, ], col = "white", 
            lwd = 2.5)
      lines(x = pred_vals, y = preds[i, ], col = sample(c("#DCBCBC", 
                                                          "#C79999", "#B97C7C", "#A25050", "#7C0000"), 
                                                        1), lwd = 2.25)
    }
  }
  else {
    polygon(c(pred_vals, rev(pred_vals)), c(cred[1, ], rev(cred[9, 
    ])), col = c_light, border = NA)
    polygon(c(pred_vals, rev(pred_vals)), c(cred[2, ], rev(cred[8, 
    ])), col = c_light_highlight, border = NA)
    polygon(c(pred_vals, rev(pred_vals)), c(cred[3, ], rev(cred[7, 
    ])), col = c_mid, border = NA)
    polygon(c(pred_vals, rev(pred_vals)), c(cred[4, ], rev(cred[6, 
    ])), col = c_mid_highlight, border = NA)
    lines(pred_vals, cred[5, ], col = c_dark, lwd = 2.5)
  }
  box(bty = "L", lwd = 2)
  if (!missing(data_test)) {
    if (class(data_train)[1] == "list") {
      data_train <- data.frame(series = data_train$series, 
                               y = data_train$y, time = data_train$time)
      data_test <- data.frame(series = data_test$series, 
                              y = data_test$y, time = data_test$time)
    }
    last_train <- (NROW(data_train)/NCOL(object$ytimes))
    if (!realisations) {
      polygon(c(pred_vals[1:(NROW(data_train)/NCOL(object$ytimes))], 
                rev(pred_vals[1:(NROW(data_train)/NCOL(object$ytimes))])), 
              c(cred[1, 1:(NROW(data_train)/NCOL(object$ytimes))], 
                rev(cred[9, 1:(NROW(data_train)/NCOL(object$ytimes))])), 
              col = "grey70", border = NA)
      lines(pred_vals[1:(NROW(data_train)/NCOL(object$ytimes))], 
            cred[5, 1:(NROW(data_train)/NCOL(object$ytimes))], 
            col = "grey70", lwd = 2.5)
    }
    points(dplyr::bind_rows(data_train, data_test) %>% dplyr::filter(series == 
                                                                       s_name) %>% dplyr::select(time, y) %>% dplyr::distinct() %>% 
             dplyr::arrange(time) %>% dplyr::pull(y), pch = 16, 
           col = "white", cex = 0.8)
    points(dplyr::bind_rows(data_train, data_test) %>% dplyr::filter(series == 
                                                                       s_name) %>% dplyr::select(time, y) %>% dplyr::distinct() %>% 
             dplyr::arrange(time) %>% dplyr::pull(y), pch = 16, 
           col = "black", cex = 0.65)
    abline(v = last_train, col = "#FFFFFF60", lwd = 2.85)
    abline(v = last_train, col = "black", lwd = 2.5, lty = "dashed")
    truth <- as.matrix(data_test %>% dplyr::filter(series == 
                                                     s_name) %>% dplyr::select(time, y) %>% dplyr::distinct() %>% 
                         dplyr::arrange(time) %>% dplyr::pull(y))
    last_train <- length(data_train %>% dplyr::filter(series == 
                                                        s_name) %>% dplyr::select(time, y) %>% dplyr::distinct() %>% 
                           dplyr::arrange(time) %>% dplyr::pull(y))
    fc <- preds[, (last_train + 1):NCOL(preds)]
    if (all(is.na(truth))) {
      score <- NULL
      message("No non-missing values in data_test$y; cannot calculate forecast score")
    }
    else {
      if (object$family %in% c("poisson", "negative binomial", 
                               "tweedie")) {
        if (max(fc, na.rm = TRUE) > 50000) {
          score <- sum(crps_mcmc_object(as.vector(truth), 
                                        fc)[, 1], na.rm = TRUE)
          message(paste0("Out of sample CRPS:\n", score))
        }
        else {
          score <- sum(drps_mcmc_object(as.vector(truth), 
                                        fc)[, 1], na.rm = TRUE)
          message(paste0("Out of sample DRPS:\n", score))
        }
      }
      else {
        score <- sum(crps_mcmc_object(as.vector(truth), 
                                      fc)[, 1], na.rm = TRUE)
        message(paste0("Out of sample CRPS:\n", score))
      }
    }
  }
  else {
    if (class(data_train)[1] == "list") {
      data_train <- data.frame(series = data_train$series, 
                               y = data_train$y, time = data_train$time)
    }
    points(data_train %>% dplyr::filter(series == s_name) %>% 
             dplyr::select(time, y) %>% dplyr::distinct() %>% 
             dplyr::arrange(time) %>% dplyr::pull(y), pch = 16, 
           col = "white", cex = 0.8)
    points(data_train %>% dplyr::filter(series == s_name) %>% 
             dplyr::select(time, y) %>% dplyr::distinct() %>% 
             dplyr::arrange(time) %>% dplyr::pull(y), pch = 16, 
           col = "black", cex = 0.65)
  }
  if (return_forecasts) {
    if (return_score) {
      if (!missing(data_test)) {
        return(list(forecast = preds[, (last_train + 
                                          1):NCOL(preds)], score = score))
      }
      else {
        return(list(forecast = preds, score = NULL))
      }
    }
    else {
      if (!missing(data_test)) {
        return(preds[, (last_train + 1):NCOL(preds)])
      }
      else {
        return(preds)
      }
    }
  }
  first_derivs <- cbind(rep(NA, NROW(preds)), t(apply(preds, 1, diff)))
  return("derivs" = first_derivs)
}
