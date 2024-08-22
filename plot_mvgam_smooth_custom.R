plot_mvgam_smooth_custom <- function (object, series = 1, smooth, residuals = FALSE, n_resid_bins = 25, 
          realisations = FALSE, n_realisations = 15, derivatives = FALSE, 
          newdata) 
{
  if (class(object) != "mvgam") {
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
  if (sign(n_resid_bins) != 1) {
    stop("argument \"n_resid_bins\" must be a positive integer", 
         call. = FALSE)
  }
  else {
    if (n_resid_bins%%1 != 0) {
      stop("argument \"n_resid_bins\" must be a positive integer", 
           call. = FALSE)
    }
  }
  s_name <- levels(object$obs_data$series)[series]
  data_train <- object$obs_data
  smooth_terms <- unlist(purrr::map(object$mgcv_model$smooth, 
                                    "label"))
  if (is.character(smooth)) {
    if (!grepl("\\(", smooth)) {
      smooth <- paste0("s(", smooth, ")")
    }
    if (!smooth %in% smooth_terms) {
      stop(smooth, " not found in smooth terms of object\nAppropriate names are: ", 
           paste(smooth_terms, collapse = ", "))
    }
    smooth_int <- which(smooth_terms == smooth)
  }
  else {
    smooth_int <- smooth
  }
  if (!object$mgcv_model$smooth[[smooth_int]]$plot.me) {
    stop(paste0("unable to plot ", object$mgcv_model$smooth[[smooth_int]]$label, 
                " (class = ", attr(object$mgcv_model$smooth[[smooth_int]], 
                                   "class")[1]), ")")
  }
  if (is.numeric(smooth)) {
    if (!smooth %in% seq_along(smooth_terms)) {
      stop(smooth, " not found in smooth terms of object")
    }
    smooth_int <- smooth
    smooth <- smooth_terms[smooth]
  }
  if (length(unlist(strsplit(smooth, ","))) > 3) {
    stop("mvgam cannot yet plot smooths of more than 3 dimensions")
  }
  smooth_labs <- do.call(rbind, lapply(seq_along(object$mgcv_model$smooth), 
                                       function(x) {
                                         data.frame(label = object$mgcv_model$smooth[[x]]$label, 
                                                    class = class(object$mgcv_model$smooth[[x]])[1])
                                       }))
  if (smooth_labs$class[smooth_int] == "random.effect") {
    stop("use function \"plot_mvgam_randomeffects\" to plot \"re\" bases")
  }
  smooth_terms <- unique(trimws(strsplit(gsub("\\+", ",", 
                                              as.character(object$mgcv_model$pred.formula)[2]), ",")[[1]]))
  smooth_terms <- smooth_terms[!grepl(",", smooth_terms)]
  smooth <- object$mgcv_model$smooth[[smooth_int]]$term
  if (length(unlist(strsplit(smooth, ","))) >= 2) {
    suppressWarnings(plot(object$mgcv_model, select = smooth_int, 
                          residuals = residuals, scheme = 2, main = "", too.far = 0, 
                          contour.col = "black", hcolors = hcl.colors(25, 
                                                                      palette = "Reds 2"), lwd = 1, seWithMean = TRUE))
    title(object$mgcv_model$smooth[[smooth_int]]$label, 
          adj = 0)
  }
  else {
    if (missing(newdata) && class(object$obs_data)[1] != 
        "list") {
      pred_dat <- data_train %>% dplyr::select(c(series, 
                                                 smooth_terms)) %>% dplyr::filter(series == s_name) %>% 
        dplyr::mutate(series = s_name)
      if (derivatives) {
        pred_dat <- pred_dat %>% dplyr::select(-smooth) %>% 
          dplyr::distinct() %>% dplyr::slice_head(n = 1) %>% 
          dplyr::bind_cols(smooth.var = seq(min(pred_dat[, 
                                                         smooth]), max(pred_dat[, smooth]), length.out = 1000))
      }
      else {
        pred_dat <- pred_dat %>% dplyr::select(-smooth) %>% 
          dplyr::distinct() %>% dplyr::slice_head(n = 1) %>% 
          dplyr::bind_cols(smooth.var = seq(min(pred_dat[, 
                                                         smooth]), max(pred_dat[, smooth]), length.out = 500))
      }
      colnames(pred_dat) <- gsub("smooth.var", smooth, 
                                 colnames(pred_dat))
    }
    else if (missing(newdata) && class(object$obs_data)[1] == 
             "list") {
      pred_dat <- vector(mode = "list")
      for (x in 1:length(data_train)) {
        if (is.matrix(data_train[[x]])) {
          pred_dat[[x]] <- matrix(0, nrow = 500, ncol = NCOL(data_train[[x]]))
        }
        else {
          pred_dat[[x]] <- rep(0, 500)
        }
      }
      names(pred_dat) <- names(object$obs_data)
      pred_dat$series <- rep((levels(data_train$series)[series]), 
                             500)
      if (!is.matrix(pred_dat[[smooth]])) {
        pred_dat[[smooth]] <- seq(min(data_train[[smooth]]), 
                                  max(data_train[[smooth]]), length.out = 500)
      }
      else {
        pred_dat[[smooth]] <- matrix(seq(min(data_train[[smooth]]), 
                                         max(data_train[[smooth]]), length.out = length(pred_dat[[smooth]])), 
                                     nrow = nrow(pred_dat[[smooth]]), ncol = ncol(pred_dat[[smooth]]))
      }
      if ("lag" %in% names(pred_dat)) {
        pred_dat[["lag"]] <- matrix(0:(NCOL(data_train$lag) - 
                                         1), nrow(pred_dat$lag), NCOL(data_train$lag), 
                                    byrow = TRUE)
      }
    }
    else {
      pred_dat <- newdata
      if (class(pred_dat)[1] != "list") {
        if (!"series" %in% colnames(pred_dat)) {
          pred_dat$series <- factor("series1")
        }
      }
      if (class(pred_dat)[1] == "list") {
        if (!"series" %in% names(pred_dat)) {
          pred_dat$series <- factor("series1")
        }
      }
    }
    suppressWarnings(Xp <- try(predict(object$mgcv_model, 
                                       newdata = pred_dat, type = "lpmatrix"), silent = TRUE))
    if (inherits(Xp, "try-error")) {
      testdat <- data.frame(series = pred_dat$series)
      terms_include <- names(object$mgcv_model$coefficients)[which(!names(object$mgcv_model$coefficients) %in% 
                                                                     "(Intercept)")]
      if (length(terms_include) > 0) {
        newnames <- vector()
        newnames[1] <- "series"
        for (i in 1:length(terms_include)) {
          testdat <- cbind(testdat, data.frame(pred_dat[[terms_include[i]]]))
          newnames[i + 1] <- terms_include[i]
        }
        colnames(testdat) <- newnames
      }
      suppressWarnings(Xp <- predict(object$mgcv_model, 
                                     newdata = testdat, type = "lpmatrix"))
    }
    keeps <- object$mgcv_model$smooth[[smooth_int]]$first.para:object$mgcv_model$smooth[[smooth_int]]$last.para
    Xp[, !seq_len(length.out = NCOL(Xp)) %in% keeps] <- 0
    if (class(pred_dat)[1] == "list") {
      if (is.matrix(pred_dat[[smooth]])) {
        pred_vals <- as.vector(as.matrix(pred_dat[[smooth]][, 
                                                            1]))
      }
      else {
        pred_vals <- as.vector(as.matrix(pred_dat[[smooth]]))
      }
    }
    else {
      pred_vals <- as.vector(as.matrix(pred_dat[, smooth]))
    }
    if (object$mgcv_model$smooth[[smooth_int]]$by != "NA") {
      by <- rep(1, length(pred_vals))
      dat <- data.frame(x = pred_vals, by = by)
      names(dat) <- c(object$mgcv_model$smooth[[smooth_int]]$term, 
                      object$mgcv_model$smooth[[smooth_int]]$by)
      Xp_term <- mgcv::PredictMat(object$mgcv_model$smooth[[smooth_int]], 
                                  dat)
      Xp[, object$mgcv_model$smooth[[smooth_int]]$first.para:object$mgcv_model$smooth[[smooth_int]]$last.para] <- Xp_term
    }
    betas <- mcmc_chains(object$model_output, "b")
    preds <- matrix(NA, nrow = NROW(betas), ncol = NROW(Xp))
    for (i in 1:NROW(betas)) {
      preds[i, ] <- (Xp %*% betas[i, ])
    }
    if (residuals) {
      suppressWarnings(Xp2 <- predict(object$mgcv_model, 
                                      newdata = object$obs_data, type = "lpmatrix"))
      if (!missing(newdata)) {
        stop("Partial residual plots not available when using newdata")
      }
      if (object$mgcv_model$smooth[[smooth_int]]$by != 
          "NA") {
        by <- rep(1, length(object$obs_data$series))
        dat <- data.frame(x = object$obs_data[[object$mgcv_model$smooth[[smooth_int]]$term]], 
                          by = by)
        names(dat) <- c(object$mgcv_model$smooth[[smooth_int]]$term, 
                        object$mgcv_model$smooth[[smooth_int]]$by)
        Xp_term <- mgcv::PredictMat(object$mgcv_model$smooth[[smooth_int]], 
                                    dat)
        Xp2[, object$mgcv_model$smooth[[smooth_int]]$first.para:object$mgcv_model$smooth[[smooth_int]]$last.para] <- Xp_term
      }
      if (class(pred_dat)[1] == "list") {
        end_train <- length(which(object$obs_data[["series"]] == 
                                    (levels(data_train$series)[series])))
      }
      else {
        end_train <- object$obs_data %>% dplyr::filter(series == 
                                                         s_name) %>% NROW()
      }
      Xp2 <- Xp2[object$ytimes[, series][1:end_train], 
      ]
      Xp2[, !grepl(paste0("(", smooth, ")"), colnames(Xp), 
                   fixed = T)] <- 0
      all_resids <- object$resids[[series]][, 1:end_train]
      partial_resids <- matrix(NA, nrow = nrow(betas), 
                               ncol = NCOL(all_resids))
      for (i in 1:NROW(betas)) {
        partial_resids[i, ] <- (Xp2 %*% betas[i, ]) + 
          all_resids[i, ]
      }
    }
    probs = c(0.05, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.95)
    cred <- sapply(1:NCOL(preds), function(n) quantile(preds[, 
                                                             n], probs = probs))
    c_light <- c("#DCBCBC")
    c_light_highlight <- c("#C79999")
    c_mid <- c("#B97C7C")
    c_mid_highlight <- c("#A25050")
    c_dark <- c("#8F2727")
    c_dark_highlight <- c("#7C0000")
    if (derivatives) {
      .pardefault <- par(no.readonly = T)
      par(.pardefault)
      par(mfrow = c(2, 1), mgp = c(2.5, 1, 0), mai = c(0.8, 
                                                       0.8, 0.4, 0.4))
      if (residuals) {
        plot(1, type = "n", bty = "L", xlab = smooth, 
             ylab = "Partial effect", xlim = c(min(pred_vals), 
                                               max(pred_vals)), ylim = c(min(min(partial_resids, 
                                                                                 min(cred) - 0.8 * sd(preds), na.rm = T)), 
                                                                         max(max(partial_resids, max(cred) + 0.8 * 
                                                                                   sd(preds), na.rm = T))))
        if (object$mgcv_model$smooth[[smooth_int]]$by != 
            "NA") {
          title(object$mgcv_model$smooth[[smooth_int]]$label, 
                adj = 0)
        }
        else {
          title(paste0("s(", smooth, ") for ", unique(pred_dat$series)), 
                adj = 0)
        }
      }
      else {
        plot(1, type = "n", bty = "L", xlab = smooth, 
             ylab = "Partial effect", xlim = c(min(pred_vals), 
                                               max(pred_vals)), ylim = c(min(cred) - 0.8 * 
                                                                           sd(preds), max(cred) + 0.8 * sd(preds)))
        if (object$mgcv_model$smooth[[smooth_int]]$by != 
            "NA") {
          title(object$mgcv_model$smooth[[smooth_int]]$label, 
                adj = 0)
        }
        else {
          title(paste0("s(", smooth, ") for ", unique(pred_dat$series)), 
                adj = 0)
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
        if (residuals) {
          sorted_x <- sort(unique(round(object$obs_data[[smooth]], 
                                        6)))
          s_name <- levels(object$obs_data$series)[series]
          obs_x <- round(data.frame(series = object$obs_data$series, 
                                    smooth_vals = object$obs_data[[smooth]]) %>% 
                           dplyr::filter(series == s_name) %>% dplyr::pull(smooth_vals), 
                         6)
          if (length(sorted_x) > n_resid_bins) {
            sorted_x <- seq(min(sorted_x), max(sorted_x), 
                            length.out = n_resid_bins)
            resid_probs <- do.call(rbind, lapply(2:n_resid_bins, 
                                                 function(i) {
                                                   quantile(as.vector(partial_resids[, 
                                                                                     which(obs_x <= sorted_x[i] & obs_x > 
                                                                                             sorted_x[i - 1])]), probs = probs)
                                                 }))
            resid_probs <- rbind(quantile(as.vector(partial_resids[, 
                                                                   which(obs_x == sorted_x[1])]), probs = probs), 
                                 resid_probs)
          }
          else {
            resid_probs <- do.call(rbind, lapply(sorted_x, 
                                                 function(i) {
                                                   quantile(as.vector(partial_resids[, 
                                                                                     which(obs_x == i)]), probs = probs)
                                                 }))
          }
          N <- length(sorted_x)
          idx <- rep(1:N, each = 2)
          repped_x <- rep(sorted_x, each = 2)
          x <- sapply(1:length(idx), function(k) if (k%%2 == 
                                                     0) 
            repped_x[k] + min(diff(sorted_x))/2
            else repped_x[k] - min(diff(sorted_x))/2)
          rect(xleft = x[seq(1, N * 2, by = 2)], xright = x[seq(2, 
                                                                N * 2, by = 2)], ytop = resid_probs[, 9], 
               ybottom = resid_probs[, 1], col = c_light, 
               border = "transparent")
          rect(xleft = x[seq(1, N * 2, by = 2)], xright = x[seq(2, 
                                                                N * 2, by = 2)], ytop = resid_probs[, 8], 
               ybottom = resid_probs[, 2], col = c_light_highlight, 
               border = "transparent")
          rect(xleft = x[seq(1, N * 2, by = 2)], xright = x[seq(2, 
                                                                N * 2, by = 2)], ytop = resid_probs[, 7], 
               ybottom = resid_probs[, 3], col = c_mid, 
               border = "transparent")
          rect(xleft = x[seq(1, N * 2, by = 2)], xright = x[seq(2, 
                                                                N * 2, by = 2)], ytop = resid_probs[, 6], 
               ybottom = resid_probs[, 4], col = c_mid_highlight, 
               border = "transparent")
          for (k in 1:N) {
            lines(x = c(x[seq(1, N * 2, by = 2)][k], 
                        x[seq(2, N * 2, by = 2)][k]), y = c(resid_probs[k, 
                                                                        5], resid_probs[k, 5]), col = c_dark, 
                  lwd = 2)
          }
          polygon(c(pred_vals, rev(pred_vals)), c(cred[1, 
          ], rev(cred[9, ])), col = rgb(red = 0, green = 0, 
                                        blue = 0, alpha = 30, maxColorValue = 200), 
          border = NA)
          lines(pred_vals, cred[5, ], col = rgb(red = 0, 
                                                green = 0, blue = 0, alpha = 45, maxColorValue = 200), 
                lwd = 3)
          box(bty = "L", lwd = 2)
        }
        else {
          polygon(c(pred_vals, rev(pred_vals)), c(cred[1, 
          ], rev(cred[9, ])), col = c_light, border = NA)
          polygon(c(pred_vals, rev(pred_vals)), c(cred[2, 
          ], rev(cred[8, ])), col = c_light_highlight, 
          border = NA)
          polygon(c(pred_vals, rev(pred_vals)), c(cred[3, 
          ], rev(cred[7, ])), col = c_mid, border = NA)
          polygon(c(pred_vals, rev(pred_vals)), c(cred[4, 
          ], rev(cred[6, ])), col = c_mid_highlight, 
          border = NA)
          lines(pred_vals, cred[5, ], col = c_dark, 
                lwd = 2.5)
        }
      }
      box(bty = "L", lwd = 2)
      if (class(object$obs_data)[1] == "list") {
        rug((as.vector(as.matrix(pred_dat[[smooth]])))[which(pred_dat[["series"]] == 
                                                               levels(pred_dat[["series"]])[series])], lwd = 1.75, 
            ticksize = 0.025, col = c_mid_highlight)
      }
      else {
        rug((as.vector(as.matrix(data_train[, smooth])))[which(data_train$series == 
                                                                 levels(data_train$series)[series])], lwd = 1.75, 
            ticksize = 0.025, col = c_mid_highlight)
      }
      first_derivs <- cbind(rep(NA, NROW(preds)), t(apply(preds, 
                                                          1, diff)))
      cred <- sapply(1:NCOL(first_derivs), function(n) quantile(first_derivs[, 
                                                                             n], probs = probs, na.rm = T))
      plot(1, type = "n", bty = "L", xlab = smooth, ylab = "1st derivative", 
           xlim = c(min(pred_vals), max(pred_vals)), ylim = c(min(cred, 
                                                                  na.rm = T) - sd(first_derivs, na.rm = T), 
                                                              max(cred, na.rm = T) + sd(first_derivs, na.rm = T)))
      if (realisations) {
        for (i in 1:n_realisations) {
          lines(x = pred_vals, y = first_derivs[i, ], 
                col = "white", lwd = 2.5)
          lines(x = pred_vals, y = first_derivs[i, ], 
                col = sample(c("#DCBCBC", "#C79999", "#B97C7C", 
                               "#A25050", "#7C0000"), 1), lwd = 2.25)
        }
      }
      else {
        polygon(c(pred_vals, rev(pred_vals)), c(cred[1, 
        ], rev(cred[9, ])), col = c_light, border = NA)
        polygon(c(pred_vals, rev(pred_vals)), c(cred[2, 
        ], rev(cred[8, ])), col = c_light_highlight, 
        border = NA)
        polygon(c(pred_vals, rev(pred_vals)), c(cred[3, 
        ], rev(cred[7, ])), col = c_mid, border = NA)
        polygon(c(pred_vals, rev(pred_vals)), c(cred[4, 
        ], rev(cred[6, ])), col = c_mid_highlight, 
        border = NA)
        lines(pred_vals, cred[5, ], col = c_dark, lwd = 2.5)
      }
      box(bty = "L", lwd = 2)
      abline(h = 0, lty = "dashed", lwd = 2)
      invisible()
      par(.pardefault)
    }
    else {
      if (residuals) {
        plot(1, type = "n", bty = "L", xlab = smooth, 
             ylab = "Partial effect", xlim = c(min(pred_vals), 
                                               max(pred_vals)), ylim = c(min(min(partial_resids, 
                                                                                 min(cred) - 0.8 * sd(preds), na.rm = T)), 
                                                                         max(max(partial_resids, max(cred) + 0.8 * 
                                                                                   sd(preds), na.rm = T))))
        if (object$mgcv_model$smooth[[smooth_int]]$by != 
            "NA") {
          title(object$mgcv_model$smooth[[smooth_int]]$label, 
                adj = 0)
        }
        else {
          title(paste0("s(", smooth, ") for ", unique(pred_dat$series)), 
                adj = 0)
        }
        sorted_x <- sort(unique(round(object$obs_data[[smooth]], 
                                      6)))
        s_name <- levels(object$obs_data$series)[series]
        obs_x <- round(data.frame(series = object$obs_data$series, 
                                  smooth_vals = object$obs_data[[smooth]]) %>% 
                         dplyr::filter(series == s_name) %>% dplyr::pull(smooth_vals), 
                       6)
        if (length(sorted_x) > n_resid_bins) {
          sorted_x <- seq(min(sorted_x), max(sorted_x), 
                          length.out = n_resid_bins)
          resid_probs <- do.call(rbind, lapply(2:n_resid_bins, 
                                               function(i) {
                                                 quantile(as.vector(partial_resids[, which(obs_x <= 
                                                                                             sorted_x[i] & obs_x > sorted_x[i - 1])]), 
                                                          probs = probs)
                                               }))
          resid_probs <- rbind(quantile(as.vector(partial_resids[, 
                                                                 which(obs_x == sorted_x[1])]), probs = probs), 
                               resid_probs)
        }
        else {
          resid_probs <- do.call(rbind, lapply(sorted_x, 
                                               function(i) {
                                                 quantile(as.vector(partial_resids[, which(obs_x == 
                                                                                             i)]), probs = probs)
                                               }))
        }
        N <- length(sorted_x)
        idx <- rep(1:N, each = 2)
        repped_x <- rep(sorted_x, each = 2)
        x <- sapply(1:length(idx), function(k) if (k%%2 == 
                                                   0) 
          repped_x[k] + min(diff(sorted_x))/2
          else repped_x[k] - min(diff(sorted_x))/2)
        rect(xleft = x[seq(1, N * 2, by = 2)], xright = x[seq(2, 
                                                              N * 2, by = 2)], ytop = resid_probs[, 9], 
             ybottom = resid_probs[, 1], col = c_light, 
             border = "transparent")
        rect(xleft = x[seq(1, N * 2, by = 2)], xright = x[seq(2, 
                                                              N * 2, by = 2)], ytop = resid_probs[, 8], 
             ybottom = resid_probs[, 2], col = c_light_highlight, 
             border = "transparent")
        rect(xleft = x[seq(1, N * 2, by = 2)], xright = x[seq(2, 
                                                              N * 2, by = 2)], ytop = resid_probs[, 7], 
             ybottom = resid_probs[, 3], col = c_mid, border = "transparent")
        rect(xleft = x[seq(1, N * 2, by = 2)], xright = x[seq(2, 
                                                              N * 2, by = 2)], ytop = resid_probs[, 6], 
             ybottom = resid_probs[, 4], col = c_mid_highlight, 
             border = "transparent")
        for (k in 1:N) {
          lines(x = c(x[seq(1, N * 2, by = 2)][k], x[seq(2, 
                                                         N * 2, by = 2)][k]), y = c(resid_probs[k, 
                                                                                                5], resid_probs[k, 5]), col = c_dark, lwd = 2)
        }
        polygon(c(pred_vals, rev(pred_vals)), c(cred[1, 
        ], rev(cred[9, ])), col = rgb(red = 0, green = 0, 
                                      blue = 0, alpha = 30, maxColorValue = 200), 
        border = NA)
        lines(pred_vals, cred[5, ], col = rgb(red = 0, 
                                              green = 0, blue = 0, alpha = 45, maxColorValue = 200), 
              lwd = 3)
        box(bty = "L", lwd = 2)
      }
      else {
        plot(1, type = "n", bty = "L", xlab = smooth, 
             ylab = "Partial effect", xlim = c(min(pred_vals), 
                                               max(pred_vals)), ylim = c(min(cred) - 0.8 * 
                                                                           sd(preds), max(cred) + 0.8 * sd(preds)))
        if (object$mgcv_model$smooth[[smooth_int]]$by != 
            "NA") {
          title(object$mgcv_model$smooth[[smooth_int]]$label, 
                adj = 0)
        }
        else {
          title(paste0("s(", smooth, ") for ", unique(pred_dat$series)), 
                adj = 0)
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
          polygon(c(pred_vals, rev(pred_vals)), c(cred[1, 
          ], rev(cred[9, ])), col = c_light, border = NA)
          polygon(c(pred_vals, rev(pred_vals)), c(cred[2, 
          ], rev(cred[8, ])), col = c_light_highlight, 
          border = NA)
          polygon(c(pred_vals, rev(pred_vals)), c(cred[3, 
          ], rev(cred[7, ])), col = c_mid, border = NA)
          polygon(c(pred_vals, rev(pred_vals)), c(cred[4, 
          ], rev(cred[6, ])), col = c_mid_highlight, 
          border = NA)
          lines(pred_vals, cred[5, ], col = c_dark, 
                lwd = 2.5)
        }
        box(bty = "L", lwd = 2)
      }
      if (class(object$obs_data)[1] == "list") {
        rug((as.vector(as.matrix(data_train[[smooth]])))[which(data_train$series == 
                                                                 levels(data_train$series)[series])], lwd = 1.75, 
            ticksize = 0.025, col = c_mid_highlight)
      }
      else {
        rug((as.vector(as.matrix(data_train[, smooth])))[which(data_train$series == 
                                                                 levels(data_train$series)[series])], lwd = 1.75, 
            ticksize = 0.025, col = c_mid_highlight)
      }
    }
  }
  
  return("derivs" = first_derivs)
}

