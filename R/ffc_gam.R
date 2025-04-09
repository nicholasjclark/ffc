#' Fit a functional time series model using dynamic functional coefficients
#'
#' Fit Generalized Additive Models that can include time-varying (dynamic)
#' functions
#'
#' @importFrom mgcv gam bam
#' @inheritParams mgcv::bam
#' @param formula A GAM formula (see \code{\link[mgcv]{formula.gam}} and also
#' \code{\link[mgcv]{gam.models}}). This is exactly like the formula for a
#' GLM except that smooth terms, [fts()], \code{\link[mgcv]{s}()} and
#' \code{\link[mgcv]{te}()} can be added to the right hand side to specify
#' that the linear predictor depends on smooth functions of predictors
#' (or linear functionals of these).
#' @param engine `character` string specifying which \pkg{mgcv} interface to use
#' for fitting the model.
#' @param time `character` specifying which variable in `data` represents the
#' the time ordering of the observations
#' @param ... 	other arguments to pass to either \code{\link[mgcv]{gam}}
#' @rdname ffc_gam
#' @details This function will update the supplied `formula` to ensure any time-varying
#' functionals (supplied through [fts()] terms in the formula right hand side) are
#' appropriately incorporated into the model. It then passes the updated model and data
#' objects to the specified `engine` for model fitting
#' @return An object of class `ffc_gam`, which inherits from objects of class `gam` or
#' `bam`
#' @seealso [fts()], \code{\link[mgcv]{gam}}, \code{\link[mgcv]{bam}}
#' @author Nicholas J Clark
#' @export
ffc_gam <- function(
    formula,
    family = gaussian(),
    data = list(),
    time,
    engine = c("gam", "bam"),
    ...) {
  engine <- match.arg(engine)
  orig_call <- match.call()
  orig_formula <- formula
  if (!exists(time, data)) {
    stop(
      paste0(
        "the variable '",
        time,
        "' cannot be found in data"
      ),
      call. = FALSE
    )
  }

  # Update formula and data by checking for any fts() terms
  interpreted <- interpret_ffc(
    formula = formula,
    data = data,
    time_var = time
  )

  # Fit the model using the specified engine
  dots <- list(...)
  fit_args <- list(
    formula = interpreted$formula,
    family = family,
    data = interpreted$data,
    na.action = "na.fail",
    ...
  )
  out <- do.call(engine, fit_args)

  # Update the object and return
  out$call <- orig_call
  out$orig_formula <- orig_formula
  out$time_var <- time
  out$fts_smooths <- interpreted$fts_smooths
  out$gam_init <- interpreted$gam_init
  out <- update_mod_data(
    gam_object = out,
    fts_smooths = interpreted$fts_smooths,
    data = data
  )

  class(out) <- c("ffc_gam", class(out))
  return(out)
}

#' Ensure all terms are included in the stored model data
#' @noRd
update_mod_data <- function(
    gam_object,
    fts_smooths,
    data) {
  # All terms included in fts smooths
  all_terms <- unique(
    unlist(purrr::map(fts_smooths, "term"))
  )

  all_terms <- unique(
    c(
      all_terms,
      unlist(purrr::map(fts_smooths, "by"))
    )
  )

  # Any offset terms
  termlabs <- attr(terms.formula(gam_object$formula, keep.order = TRUE), "term.labels")

  # Check for offsets as well
  off_names <- grep(
    "offset",
    rownames(attr(terms.formula(gam_object$formula), "factors")),
    value = TRUE
  )
  if (length(off_names) > 0L) {
    all_terms <- c(all_terms, strip_offset(off_names))
  }

  # Add any terms that aren't already in the model slot
  vars_to_add <- setdiff(
    all_terms,
    c("NA", colnames(gam_object$model))
  )

  if (length(vars_to_add)) {
    orig_names <- colnames(gam_object$model)
    for (i in 1:length(vars_to_add)) {
      gam_object$model <- cbind(
        gam_object$model,
        data[[vars_to_add[i]]]
      )
    }
    colnames(gam_object$model) <- c(orig_names, vars_to_add)
  }
  return(gam_object)
}

#' Strip offset names down to the actual terms
#' @noRd
strip_offset <- function(x) {
  for (i in 1:length(x)) {
    if (substr(x[i], 1, 7) == "offset(") {
      x[i] <- substr(x[i], 8, nchar(x[i]) - 1)
    }

    if (substr(x[i], 1, 4) == "log(") {
      x[i] <- substr(x[i], 5, nchar(x[i]) - 1)
    }

    if (substr(x[i], 1, 4) == "exp(") {
      x[i] <- substr(x[i], 5, nchar(x[i]) - 1)
    }

    if (substr(x[i], 1, 5) == "log2(") {
      x[i] <- substr(x[i], 6, nchar(x[i]) - 1)
    }

    if (substr(x[i], 1, 5) == "sqrt(") {
      x[i] <- substr(x[i], 6, nchar(x[i]) - 1)
    }

    if (substr(x[i], 1, 6) == "log10(") {
      x[i] <- substr(x[i], 7, nchar(x[i]) - 1)
    }
  }
  x
}

#' Generic GAM setup function
#' @importFrom stats na.fail coef gaussian model.frame quantile reformulate rnorm sd vcov
#' @importFrom utils data
#' @noRd
ffc_gam_setup <- function(
    formula,
    knots,
    family = gaussian(),
    dat = list()) {
  if (missing(knots)) {
    out <- init_gam(
      formula(formula),
      data = dat,
      family = family
    )
    attr(out, "knots") <- NULL
  } else {
    if (!is.list(knots)) {
      stop('all "knot" arguments must be supplied as lists', call. = FALSE)
    }
    out <- init_gam(
      formula(formula),
      data = dat,
      family = family,
      knots = knots
    )
    attr(out, "knots") <- knots
  }
  out
}

#' The below functions are mostly perfect copies of functions
#' written originally by Prof Simon Wood
#' All credit goes to Prof Wood and the mgcv development team.
#' They only exist in ffc_gam because of CRAN restrictions on
#' calling internal functions from other packages

#' Multivariate normal random number generator
#'
#' This function is derived from \code{mgcv::rmvn}
#'
#' @author Simon N Wood with modifications by Nicholas Clark
#' @noRd
rmvn <- function(n, mu, sig) {
  L <- mgcv::mroot(sig)
  m <- ncol(L)
  t(mu + L %*% matrix(rnorm(m * n), m, n))
}

#' Initiate a gam object without actually doing any estimation
#'
#' This function is derived from several internal mgcv functions
#'
#' @author Simon N Wood with modifications by Nicholas Clark
#' @noRd
init_gam <- function(
    formula,
    family = gaussian(),
    data = list(),
    na.action = na.omit,
    knots = NULL,
    drop.unused.levels = TRUE,
    control = mgcv::gam.control(),
    centred = TRUE,
    diagonalize = FALSE,
    sp = NULL
) {
  if (is.character(family)) family <- eval(parse(text = family))
  if (is.function(family)) family <- family()
  if (is.null(family$family)) stop("family not recognized")
  gp <- mgcv::interpret.gam(formula) # interpret the formula
  cl <- match.call() # call needed in gam object for update to work
  mf <- match.call(expand.dots = FALSE)
  mf$formula <- gp$fake.formula
  mf$family <- mf$knots <- mf$sp <- mf$file <- mf$control <-
    mf$centred <- mf$sp.prior <- mf$diagonalize <- NULL
  mf$drop.unused.levels <- drop.unused.levels
  mf[[1]] <- quote(stats::model.frame) ##as.name("model.frame")
  pmf <- mf

  pmf$formula <- gp$pf
  pmf <- eval(pmf, parent.frame())
  pterms <- attr(pmf, "terms")
  rm(pmf)

  mf <- eval(mf, parent.frame())
  if (nrow(mf) < 2) stop("Not enough (non-NA) data to do anything meaningful")
  terms <- attr(mf, "terms")

  ## summarize the *raw* input variables
  ## note can't use get_all_vars here -- buggy with matrices
  vars <- all.vars(gp$fake.formula[-2]) ## drop response here
  inp <- parse(text = paste("list(", paste(vars, collapse = ","), ")"))
  if (!is.list(data) && !is.data.frame(data)) data <- as.data.frame(data)

  dl <- eval(inp, data, parent.frame())
  if (!control$keepData) {
    rm(data)
  } ## save space
  names(dl) <- vars ## list of all variables needed
  var.summary <- variable_summary(gp$pf, dl, nrow(mf)) ## summarize the input data
  rm(dl)

  G <- gam_setup(
    gp,
    pterms = pterms,
    data = mf,
    knots = knots,
    sp = sp,
    H = NULL,
    absorb.cons = centred,
    sparse.cons = FALSE,
    select = TRUE,
    idLinksBases = TRUE,
    scale.penalty = control$scalePenalty,
    diagonal.penalty = diagonalize
  )
  G$model <- mf
  G$terms <- terms
  G$family <- family
  G$call <- cl
  G$var.summary <- var.summary

  lambda <- initial_spg(
    G$X,
    G$y,
    G$w,
    family,
    G$S,
    G$rank,
    G$off,
    offset = G$offset,
    L = G$L
  )
  jags.ini <- list()
  lam <- if (is.null(G$L)) lambda else G$L %*% lambda
  #jin <- mgcv:::jini(G,lam)
  G$formula <- formula
  G$coefficients <- rep(0, length(G$term.names))
  names(G$coefficients) <- G$term.names
  G$residuals <- rnorm(NROW(G$X))
  G$edf <- rep(1, length(G$coefficients))
  names(G$edf) <- G$term.names
  G$edf1 <- rep(1, length(G$coefficients))
  names(G$edf1) <- G$term.names
  G$sig2 <- 1
  G$rank <- ncol(G$X)
  G$Vp <- G$Ve <- diag(rep(1, length(G$coefficients)))
  G$sp <- exp(G$sp)
  G$scale.estimated <- FALSE
  G$method <- 'UBRE'
  G$pred.formula <- gp$pred.formula
  class(G) <- c('gam', 'glm', 'lm')
  G$R <- model.matrix(G)
  return(G)
}

#' Set up a gam object without actually doing any estimation
#'
#' This function is derived from \code{mgcv:::gam.setup}
#'
#' @author Simon N Wood with modifications by Nicholas Clark
#' @importFrom mgcv gam.side smoothCon get.var Rrank interpret.gam initial.sp
#' @importFrom stats .getXlevels model.matrix model.offset na.omit
#' @importFrom methods cbind2
#' @noRd
gam_setup <- function(
    formula,
    pterms,
    data = stop("No data supplied to gam_setup"),
    knots = NULL,
    sp = NULL,
    min.sp = NULL,
    H = NULL,
    absorb.cons = TRUE,
    sparse.cons = 0,
    select = FALSE,
    idLinksBases = TRUE,
    scale.penalty = TRUE,
    paraPen = NULL,
    gamm.call = FALSE,
    drop.intercept = FALSE,
    diagonal.penalty = FALSE,
    apply.by = TRUE,
    list.call = FALSE,
    modCon = 0) {
  if (inherits(formula, "split.gam.formula")) {
    split <- formula
  } else if (
    inherits(formula, "formula")
  ) {
    split <- mgcv::interpret.gam(formula)
  } else {
    stop("First argument is no sort of formula!")
  }
  if (length(split$smooth.spec) == 0) {
    if (split$pfok == 0) stop("You've got no model....")
    m <- 0
  } else {
    m <- length(split$smooth.spec)
  }
  G <- list(
    m = m,
    min.sp = min.sp,
    H = H,
    pearson.extra = 0,
    dev.extra = 0,
    n.true = -1,
    pterms = pterms
  )
  if (is.null(attr(data, "terms"))) {
    mf <- model.frame(split$pf, data, drop.unused.levels = FALSE)
  } else {
    mf <- data
  }
  G$intercept <- attr(attr(mf, "terms"), "intercept") > 0
  if (list.call) {
    offi <- attr(pterms, "offset")
    if (!is.null(offi)) {
      G$offset <- mf[[names(attr(pterms, "dataClasses"))[offi]]]
    }
  } else {
    G$offset <- model.offset(mf)
  }
  if (!is.null(G$offset)) G$offset <- as.numeric(G$offset)
  if (drop.intercept) attr(pterms, "intercept") <- 1
  X <- model.matrix(pterms, mf)
  if (drop.intercept) {
    xat <- attributes(X)
    ind <- xat$assign > 0
    X <- X[, ind, drop = FALSE]
    xat$assign <- xat$assign[ind]
    xat$dimnames[[2]] <- xat$dimnames[[2]][ind]
    xat$dim[2] <- xat$dim[2] - 1
    attributes(X) <- xat
    G$intercept <- FALSE
  }
  rownames(X) <- NULL
  G$nsdf <- ncol(X)
  G$contrasts <- attr(X, "contrasts")
  G$xlevels <- .getXlevels(pterms, mf)
  G$assign <- attr(X, "assign")
  PP <- parametric_penalty(pterms, G$assign, paraPen, sp)
  if (!is.null(PP)) {
    ind <- 1:length(PP$sp)
    if (!is.null(sp)) sp <- sp[-ind]
    if (!is.null(min.sp)) {
      PP$min.sp <- min.sp[ind]
      min.sp <- min.sp[-ind]
    }
  }
  G$smooth <- list()
  G$S <- list()
  if (gamm.call) {
    if (m > 0) for (i in 1:m) attr(split$smooth.spec[[i]], "gamm") <- TRUE
  }
  if (m > 0 && idLinksBases) {
    id.list <- list()
    for (i in 1:m) {
      if (!is.null(split$smooth.spec[[i]]$id)) {
        id <- as.character(split$smooth.spec[[i]]$id)
        if (length(id.list) && id %in% names(id.list)) {
          ni <- length(id.list[[id]]$sm.i)
          id.list[[id]]$sm.i[ni + 1] <- i
          base.i <- id.list[[id]]$sm.i[1]
          split$smooth.spec[[i]] <- clone_smooth_spec(
            split$smooth.spec[[base.i]],
            split$smooth.spec[[i]]
          )
          temp.term <- split$smooth.spec[[i]]$term
          for (j in 1:length(temp.term)) {
            id.list[[id]]$data[[j]] <- cbind(
              id.list[[id]]$data[[j]],
              mgcv::get.var(temp.term[j], data, vecMat = FALSE)
            )
          }
        } else {
          id.list[[id]] <- list(sm.i = i)
          id.list[[id]]$data <- list()
          term <- split$smooth.spec[[i]]$term
          for (j in 1:length(term)) {
            id.list[[id]]$data[[j]] <- mgcv::get.var(
              term[j],
              data,
              vecMat = FALSE
            )
          }
        }
      }
    }
  }
  G$off <- array(0, 0)
  first.para <- G$nsdf + 1
  sm <- list()
  newm <- 0
  if (m > 0) {
    for (i in 1:m) {
      id <- split$smooth.spec[[i]]$id
      if (is.null(id) || !idLinksBases) {
        sml <- mgcv::smoothCon(
          split$smooth.spec[[i]],
          data,
          knots,
          absorb.cons,
          scale.penalty = scale.penalty,
          null.space.penalty = select,
          sparse.cons = sparse.cons,
          diagonal.penalty = diagonal.penalty,
          apply.by = apply.by,
          modCon = modCon
        )
      } else {
        names(id.list[[id]]$data) <- split$smooth.spec[[i]]$term
        sml <- mgcv::smoothCon(
          split$smooth.spec[[i]],
          id.list[[id]]$data,
          knots,
          absorb.cons,
          n = nrow(data),
          dataX = data,
          scale.penalty = scale.penalty,
          null.space.penalty = select,
          sparse.cons = sparse.cons,
          diagonal.penalty = diagonal.penalty,
          apply.by = apply.by,
          modCon = modCon
        )
      }
      ind <- 1:length(sml)
      sm[ind + newm] <- sml[ind]
      newm <- newm + length(sml)
    }
  }
  G$m <- m <- newm
  if (m > 0) {
    sm <- mgcv::gam.side(sm, X, tol = .Machine$double.eps^0.5)
    if (!apply.by) {
      for (i in 1:length(sm)) {
        if (!is.null(sm[[i]]$X0)) {
          ind <- attr(sm[[i]], "del.index")
          sm[[i]]$X <- if (is.null(ind)) {
            sm[[i]]$X0
          } else {
            sm[[i]]$X0[, -ind, drop = FALSE]
          }
        }
      }
    }
  }
  idx <- list()
  L <- matrix(0, 0, 0)
  lsp.names <- sp.names <- rep("", 0)
  if (m > 0) {
    for (i in 1:m) {
      id <- sm[[i]]$id
      length.S <- if (is.null(sm[[i]]$updateS)) {
        length(sm[[i]]$S)
      } else {
        sm[[i]]$n.sp
      }
      Li <- if (is.null(sm[[i]]$L)) diag(length.S) else sm[[i]]$L
      if (length.S > 0) {
        if (length.S == 1) {
          lspn <- sm[[i]]$label
        } else {
          Sname <- names(sm[[i]]$S)
          lspn <- if (is.null(Sname)) {
            paste(sm[[i]]$label, 1:length.S, sep = "")
          } else {
            paste(sm[[i]]$label, Sname, sep = "")
          }
        }
        spn <- lspn[1:ncol(Li)]
      }
      if (is.null(id) || is.null(idx[[id]])) {
        if (!is.null(id)) {
          idx[[id]]$c <- ncol(L) + 1
          idx[[id]]$nc <- ncol(Li)
        }
        L <- rbind(
          cbind(L, matrix(0, nrow(L), ncol(Li))),
          cbind(matrix(0, nrow(Li), ncol(L)), Li)
        )
        if (length.S > 0) {
          sp.names <- c(sp.names, spn)
          lsp.names <- c(lsp.names, lspn)
        }
      } else {
        L0 <- matrix(0, nrow(Li), ncol(L))
        if (ncol(Li) > idx[[id]]$nc) {
          stop(
            "Later terms sharing an `id' can not have more smoothing parameters than the first such term"
          )
        }
        L0[, idx[[id]]$c:(idx[[id]]$c + ncol(Li) - 1)] <- Li
        L <- rbind(L, L0)
        if (length.S > 0) {
          lsp.names <- c(lsp.names, lspn)
        }
      }
    }
  }
  Xp <- NULL
  if (m > 0) {
    for (i in 1:m) {
      n.para <- ncol(sm[[i]]$X)
      sm[[i]]$first.para <- first.para
      first.para <- first.para + n.para
      sm[[i]]$last.para <- first.para - 1
      Xoff <- attr(sm[[i]]$X, "offset")
      if (!is.null(Xoff)) {
        if (is.null(G$offset)) G$offset <- Xoff else G$offset <- G$offset + Xoff
      }
      if (is.null(sm[[i]]$Xp)) {
        if (!is.null(Xp)) Xp <- cbind2(Xp, sm[[i]]$X)
      } else {
        if (is.null(Xp)) Xp <- X
        Xp <- cbind2(Xp, sm[[i]]$Xp)
        sm[[i]]$Xp <- NULL
      }
      X <- cbind2(X, sm[[i]]$X)
      sm[[i]]$X <- NULL
      G$smooth[[i]] <- sm[[i]]
    }
  }
  if (is.null(Xp)) {
    G$cmX <- colMeans(X)
  } else {
    G$cmX <- colMeans(Xp)
    qrx <- qr(Xp, LAPACK = TRUE)
    R <- qr.R(qrx)
    p <- ncol(R)
    rank <- mgcv::Rrank(R)
    QtX <- qr.qty(qrx, X)[1:rank, ]
    if (rank < p) {
      R <- R[1:rank, ]
      qrr <- qr(t(R), tol = 0)
      R <- qr.R(qrr)
      G$P <- forwardsolve(t(R), QtX)
    } else {
      G$P <- backsolve(R, QtX)
    }
    if (rank < p) {
      G$P <- qr.qy(qrr, rbind(G$P, matrix(0, p - rank, p)))
    }
    G$P[qrx$pivot, ] <- G$P
  }
  G$X <- X
  rm(X)
  n.p <- ncol(G$X)
  if (!is.null(sp)) {
    ok <- TRUE
    if (length(sp) < ncol(L)) {
      warning("Supplied smoothing parameter vector is too short - ignored.")
      ok <- FALSE
    }
    if (sum(is.na(sp))) {
      warning("NA's in supplied smoothing parameter vector - ignoring.")
      ok <- FALSE
    }
  } else {
    ok <- FALSE
  }
  G$sp <- if (ok) sp[1:ncol(L)] else rep(-1, ncol(L))
  names(G$sp) <- sp.names
  k <- 1
  if (m > 0) {
    for (i in 1:m) {
      id <- sm[[i]]$id
      if (is.null(sm[[i]]$L)) Li <- diag(length(sm[[i]]$S)) else Li <- sm[[i]]$L
      if (is.null(id)) {
        spi <- sm[[i]]$sp
        if (!is.null(spi)) {
          if (length(spi) != ncol(Li)) {
            stop(
              "incorrect number of smoothing parameters supplied for a smooth term"
            )
          }
          G$sp[k:(k + ncol(Li) - 1)] <- spi
        }
        k <- k + ncol(Li)
      } else {
        spi <- sm[[i]]$sp
        if (is.null(idx[[id]]$sp.done)) {
          if (!is.null(spi)) {
            if (length(spi) != ncol(Li)) {
              stop(
                "incorrect number of smoothing parameters supplied for a smooth term"
              )
            }
            G$sp[idx[[id]]$c:(idx[[id]]$c + idx[[id]]$nc - 1)] <- spi
          }
          idx[[id]]$sp.done <- TRUE
          k <- k + idx[[id]]$nc
        }
      }
    }
  }
  k <- 1
  if (length(idx)) for (i in 1:length(idx)) idx[[i]]$sp.done <- FALSE
  if (m > 0) {
    for (i in 1:m) {
      id <- sm[[i]]$id
      if (!is.null(id)) {
        if (idx[[id]]$nc > 0) {
          G$smooth[[i]]$sp <- G$sp[
            idx[[id]]$c:(idx[[id]]$c +
              idx[[id]]$nc -
              1)
          ]
        }
        if (!idx[[id]]$sp.done) {
          idx[[id]]$sp.done <- TRUE
          k <- k + idx[[id]]$nc
        }
      } else {
        if (is.null(sm[[i]]$L)) {
          nc <- length(sm[[i]]$S)
        } else {
          nc <- ncol(sm[[i]]$L)
        }
        if (nc > 0) G$smooth[[i]]$sp <- G$sp[k:(k + nc - 1)]
        k <- k + nc
      }
    }
  }
  if (!is.null(min.sp)) {
    if (length(min.sp) < nrow(L)) stop("length of min.sp is wrong.")
    if (nrow(L) > 0) min.sp <- min.sp[1:nrow(L)]
    if (sum(is.na(min.sp))) stop("NA's in min.sp.")
    if (sum(min.sp < 0)) stop("elements of min.sp must be non negative.")
  }
  k.sp <- 0
  G$rank <- array(0, 0)
  if (m > 0) {
    for (i in 1:m) {
      sm <- G$smooth[[i]]
      if (length(sm$S) > 0) {
        for (j in 1:length(sm$S)) {
          k.sp <- k.sp + 1
          G$off[k.sp] <- sm$first.para
          G$S[[k.sp]] <- sm$S[[j]]
          G$rank[k.sp] <- sm$rank[j]
          if (!is.null(min.sp)) {
            if (is.null(H)) H <- matrix(0, n.p, n.p)
            H[sm$first.para:sm$last.para, sm$first.para:sm$last.para] <- H[
              sm$first.para:sm$last.para,
              sm$first.para:sm$last.para
            ] +
              min.sp[k.sp] *
                sm$S[[j]]
          }
        }
      }
    }
  }
  if (!is.null(PP)) {
    L <- rbind(
      cbind(L, matrix(0, nrow(L), ncol(PP$L))),
      cbind(matrix(0, nrow(PP$L), ncol(L)), PP$L)
    )
    G$off <- c(PP$off, G$off)
    G$S <- c(PP$S, G$S)
    G$rank <- c(PP$rank, G$rank)
    G$sp <- c(PP$sp, G$sp)
    lsp.names <- c(PP$full.sp.names, lsp.names)
    G$n.paraPen <- length(PP$off)
    if (!is.null(PP$min.sp)) {
      if (is.null(H)) H <- matrix(0, n.p, n.p)
      for (i in 1:length(PP$S)) {
        ind <- PP$off[i]:(PP$off[i] + ncol(PP$S[[i]]) - 1)
        H[ind, ind] <- H[ind, ind] + PP$min.sp[i] * PP$S[[i]]
      }
    }
  } else {
    G$n.paraPen <- 0
  }
  fix.ind <- G$sp >= 0
  if (sum(fix.ind)) {
    lsp0 <- G$sp[fix.ind]
    ind <- lsp0 == 0
    ef0 <- indi <- (1:length(ind))[ind]
    if (length(indi) > 0) {
      for (i in 1:length(indi)) {
        ii <- G$off[i]:(G$off[i] + ncol(G$S[[i]]) - 1)
        ef0[i] <- norm(G$X[, ii], type = "F")^2 /
          norm(G$S[[i]], type = "F") *
          .Machine$double.eps *
          0.1
      }
    }
    lsp0[!ind] <- log(lsp0[!ind])
    lsp0[ind] <- log(ef0)
    lsp0 <- as.numeric(L[, fix.ind, drop = FALSE] %*% lsp0)
    L <- L[, !fix.ind, drop = FALSE]
    G$sp <- G$sp[!fix.ind]
  } else {
    lsp0 <- rep(0, nrow(L))
  }
  G$H <- H
  if (ncol(L) == nrow(L) && !sum(L != diag(ncol(L)))) L <- NULL
  G$L <- L
  G$lsp0 <- lsp0
  names(G$lsp0) <- lsp.names
  if (absorb.cons == FALSE) {
    G$C <- matrix(0, 0, n.p)
    if (m > 0) {
      for (i in 1:m) {
        if (is.null(G$smooth[[i]]$C)) {
          n.con <- 0
        } else {
          n.con <- nrow(G$smooth[[i]]$C)
        }
        C <- matrix(0, n.con, n.p)
        C[, G$smooth[[i]]$first.para:G$smooth[[i]]$last.para] <- G$smooth[[i]]$C
        G$C <- rbind(G$C, C)
        G$smooth[[i]]$C <- NULL
      }
      rm(C)
    }
  }
  G$y <- drop(data[[split$response]])
  ydim <- dim(G$y)
  if (!is.null(ydim) && length(ydim) < 2) dim(G$y) <- NULL
  G$n <- nrow(data)
  if (is.null(data$"(weights)")) G$w <- rep(1, G$n) else G$w <- data$"(weights)"
  if (G$nsdf > 0) {
    term.names <- colnames(G$X)[1:G$nsdf]
  } else {
    term.names <- array("", 0)
  }
  n.smooth <- length(G$smooth)
  n.sp0 <- 0
  if (n.smooth) {
    for (i in 1:n.smooth) {
      k <- 1
      jj <- G$smooth[[i]]$first.para:G$smooth[[i]]$last.para
      if (G$smooth[[i]]$df > 0) {
        for (j in jj) {
          term.names[j] <- paste(
            G$smooth[[i]]$label,
            ".",
            as.character(k),
            sep = ""
          )
          k <- k + 1
        }
      }
      n.sp <- length(G$smooth[[i]]$S)
      if (n.sp) {
        G$smooth[[i]]$first.sp <- n.sp0 + 1
        n.sp0 <- G$smooth[[i]]$last.sp <- n.sp0 + n.sp
      }
      if (!is.null(G$smooth[[i]]$g.index)) {
        if (is.null(G$g.index)) G$g.index <- rep(FALSE, n.p)
        G$g.index[jj] <- G$smooth[[i]]$g.index
      }
    }
  }
  G$term.names <- term.names
  G$pP <- PP
  G
}

#' Initiate parametric term penalties
#'
#' This function is derived from \code{mgcv:::parametricPenalty}
#'
#' @author Simon N Wood with modifications by Nicholas Clark
#' @noRd
parametric_penalty <- function(pterms, assign, paraPen, sp0) {
  S <- list()
  off <- rep(0, 0)
  rank <- rep(0, 0)
  sp <- rep(0, 0)
  full.sp.names <- rep("", 0)
  L <- matrix(0, 0, 0)
  k <- 0
  tind <- unique(assign)
  n.t <- length(tind)
  if (n.t > 0) {
    for (j in 1:n.t) {
      if (tind[j] > 0) {
        term.label <- attr(pterms[tind[j]], "term.label")
        P <- paraPen[[term.label]]
        if (!is.null(P)) {
          ind <- (1:length(assign))[assign == tind[j]]
          Li <- P$L
          P$L <- NULL
          spi <- P$sp
          P$sp <- NULL
          ranki <- P$rank
          P$rank <- NULL
          np <- length(P)
          if (!is.null(ranki) && length(ranki) != np) {
            stop("`rank' has wrong length in `paraPen'")
          }
          if (np) {
            for (i in 1:np) {
              k <- k + 1
              S[[k]] <- P[[i]]
              off[k] <- min(ind)
              if (ncol(P[[i]]) != nrow(P[[i]]) || nrow(P[[i]]) != length(ind)) {
                stop(" a parametric penalty has wrong dimension")
              }
              if (is.null(ranki)) {
                ev <- eigen(S[[k]], symmetric = TRUE, only.values = TRUE)$values
                rank[k] <- sum(ev > max(ev) * .Machine$double.eps * 10)
              } else {
                rank[k] <- ranki[i]
              }
            }
          }
          if (np) {
            if (is.null(Li)) Li <- diag(np)
            if (nrow(Li) != np) stop("L has wrong dimension in `paraPen'")
            L <- rbind(
              cbind(L, matrix(0, nrow(L), ncol(Li))),
              cbind(matrix(0, nrow(Li), ncol(L)), Li)
            )
            ind <- (length(sp) + 1):(length(sp) + ncol(Li))
            ind2 <- (length(sp) + 1):(length(sp) + nrow(Li))
            if (is.null(spi)) {
              sp[ind] <- -1
            } else {
              if (length(spi) != ncol(Li)) {
                stop("`sp' dimension wrong in `paraPen'")
              }
              sp[ind] <- spi
            }
            if (length(ind) > 1) {
              names(sp)[ind] <- paste(
                term.label,
                ind -
                  ind[1] +
                  1,
                sep = ""
              )
            } else {
              names(sp)[ind] <- term.label
            }
            if (length(ind2) > 1) {
              full.sp.names[ind2] <- paste(
                term.label,
                ind2 - ind2[1] + 1,
                sep = ""
              )
            } else {
              full.sp.names[ind2] <- term.label
            }
          }
        }
      }
    }
  }
  if (k == 0) {
    return(NULL)
  }
  if (!is.null(sp0)) {
    if (length(sp0) < length(sp)) stop("`sp' too short")
    sp0 <- sp0[1:length(sp)]
    sp[sp < 0] <- sp0[sp < 0]
  }
  list(
    S = S,
    off = off,
    sp = sp,
    L = L,
    rank = rank,
    full.sp.names = full.sp.names
  )
}

#' Clone a smooth specification
#'
#' This function is derived from \code{mgcv:::clone.smooth.spec}
#'
#' @author Simon N Wood with modifications by Nicholas Clark
#' @noRd
clone_smooth_spec <- function(specb, spec) {
  if (specb$dim != spec$dim) {
    stop("`id' linked smooths must have same number of arguments")
  }
  if (inherits(specb, c("tensor.smooth.spec", "t2.smooth.spec"))) {
    specb$term <- spec$term
    specb$label <- spec$label
    specb$by <- spec$by
    k <- 1
    for (i in 1:length(specb$margin)) {
      if (is.null(spec$margin)) {
        for (j in 1:length(specb$margin[[i]]$term)) {
          specb$margin[[i]]$term[j] <- spec$term[k]
          k <- k + 1
        }
        specb$margin[[i]]$label <- ""
      } else {
        specb$margin[[i]]$term <- spec$margin[[i]]$term
        specb$margin[[i]]$label <- spec$margin[[i]]$label
        specb$margin[[i]]$xt <- spec$margin[[i]]$xt
      }
    }
  } else {
    specb$term <- spec$term
    specb$label <- spec$label
    specb$by <- spec$by
    specb$xt <- spec$xt
  }
  specb
}

#' Summarize all the variables in a list of variables
#'
#' This function is derived from \code{mgcv:::variable.summary}
#'
#' @author Simon N Wood with modifications by Nicholas Clark
#' @noRd
variable_summary <- function(pf, dl, n) {
  v.n <- length(dl)
  v.name <- v.name1 <- names(dl)
  if (v.n) {
    k <- 0
    for (i in 1:v.n) {
      if (length(dl[[i]]) >= n) {
        k <- k + 1
        v.name[k] <- v.name1[i]
      }
    }
    if (k > 0) v.name <- v.name[1:k] else v.name <- rep("", k)
  }
  p.name <- all.vars(pf[-2])
  vs <- list()
  v.n <- length(v.name)
  if (v.n > 0) {
    for (i in 1:v.n) {
      if (v.name[i] %in% p.name) para <- TRUE else para <- FALSE
      if (para && is.matrix(dl[[v.name[i]]]) && ncol(dl[[v.name[i]]]) > 1) {
        x <- matrix(
          apply(
            dl[[v.name[i]]],
            2,
            quantile,
            probs = 0.5,
            type = 3,
            na.rm = TRUE
          ),
          1,
          ncol(dl[[v.name[i]]])
        )
      } else {
        x <- dl[[v.name[i]]]
        if (is.character(x)) x <- as.factor(x)
        if (is.factor(x)) {
          x <- x[!is.na(x)]
          lx <- levels(x)
          freq <- tabulate(x)
          ii <- min((1:length(lx))[freq == max(freq)])
          x <- factor(lx[ii], levels = lx)
        } else {
          x <- as.numeric(x)
          x <- c(
            min(x, na.rm = TRUE),
            as.numeric(quantile(x, probs = 0.5, type = 3, na.rm = TRUE)),
            max(x, na.rm = TRUE)
          )
        }
      }
      vs[[v.name[i]]] <- x
    }
  }
  vs
}

#' Initialize smoothing parameters
#'
#' This function is derived from \code{mgcv:::initial.spg}
#'
#' @author Simon N Wood with modifications by Nicholas Clark
#' @importFrom stats lm
#' @noRd
initial_spg <- function(
    x,
    y,
    weights,
    family,
    S,
    rank,
    off,
    offset = NULL,
    L = NULL,
    lsp0 = NULL,
    type = 1,
    start = NULL,
    mustart = NULL,
    etastart = NULL,
    E = NULL,
    ...) {
  if (length(S) == 0) {
    return(rep(0, 0))
  }
  nobs <- nrow(x)
  if (is.null(mustart)) mukeep <- NULL else mukeep <- mustart
  eval(family$initialize)
  if (inherits(family, "general.family")) {
    lbb <- family$ll(
      y,
      x,
      start,
      weights,
      family,
      offset = offset,
      deriv = 1
    )$lbb
    pcount <- rep(0, ncol(lbb))
    for (i in 1:length(S)) {
      ind <- off[i]:(off[i] + ncol(S[[i]]) - 1)
      dlb <- -diag(lbb[ind, ind, drop = FALSE])
      indp <- rowSums(abs(S[[i]])) > max(S[[i]]) * .Machine$double.eps^0.75 &
        dlb != 0
      ind <- ind[indp]
      pcount[ind] <- pcount[ind] + 1
    }
    lambda <- rep(0, length(S))
    for (i in 1:length(S)) {
      ind <- off[i]:(off[i] + ncol(S[[i]]) - 1)
      lami <- 1
      dlb <- abs(diag(lbb[ind, ind, drop = FALSE]))
      dS <- diag(S[[i]])
      pc <- pcount[ind]
      ind <- rowSums(abs(S[[i]])) > max(S[[i]]) * .Machine$double.eps^0.75 &
        dlb != 0
      dlb <- dlb[ind] / pc[ind]
      dS <- dS[ind]
      rm <- max(length(dS) / rank[i], 1)
      while (
        sqrt(
          mean(dlb / (dlb + lami * dS * rm)) *
            mean(dlb) /
            mean(
              dlb +
                lami * dS * rm
            )
        ) >
          0.4
      ) {
        lami <- lami * 5
      }
      while (
        sqrt(
          mean(dlb / (dlb + lami * dS * rm)) *
            mean(dlb) /
            mean(
              dlb +
                lami * dS * rm
            )
        ) <
          0.4
      ) {
        lami <- lami / 5
      }
      lambda[i] <- lami
    }
  } else {
    if (is.null(mukeep)) {
      if (!is.null(start)) etastart <- drop(x %*% start)
      if (!is.null(etastart)) mustart <- family$linkinv(etastart)
    } else {
      mustart <- mukeep
    }
    if (inherits(family, "extended.family")) {
      theta <- family$getTheta()
      Ddo <- family$Dd(y, mustart, theta, weights)
      mu.eta2 <- family$mu.eta(family$linkfun(mustart))^2
      w <- 0.5 * as.numeric(Ddo$Dmu2 * mu.eta2)
      if (any(w < 0)) w <- 0.5 * as.numeric(Ddo$EDmu2 * mu.eta2)
    } else {
      w <- as.numeric(
        weights *
          family$mu.eta(family$linkfun(mustart))^2 /
          family$variance(mustart)
      )
    }
    w <- sqrt(w)
    if (type == 1) {
      lambda <- mgcv::initial.sp(w * x, S, off)
    } else {
      csX <- colSums((w * x)^2)
      lambda <- rep(0, length(S))
      for (i in 1:length(S)) {
        ind <- off[i]:(off[i] + ncol(S[[i]]) - 1)
        lambda[i] <- sum(csX[ind]) / sqrt(sum(S[[i]]^2))
      }
    }
  }
  if (!is.null(L)) {
    lsp <- log(lambda)
    if (is.null(lsp0)) lsp0 <- rep(0, nrow(L))
    lsp <- as.numeric(coef(lm(lsp ~ L - 1 + offset(lsp0))))
    lambda <- exp(lsp)
  }
  lambda
}

#' Initialize a gam object using a list of formulae
#'
#' This function is derived from \code{mgcv:::gam.setup.list}
#'
#' @author Simon N Wood with modifications by Nicholas Clark
#' @noRd
gam_setup.list <- function(
    formula,
    pterms,
    data = stop("No data supplied to gam.setup"),
    knots = NULL,
    sp = NULL,
    min.sp = NULL,
    H = NULL,
    absorb.cons = TRUE,
    sparse.cons = 0,
    select = FALSE,
    idLinksBases = TRUE,
    scale.penalty = TRUE,
    paraPen = NULL,
    gamm.call = FALSE,
    drop.intercept = NULL,
    apply.by = TRUE,
    modCon = 0) {
  d <- length(pterms)
  if (is.null(drop.intercept)) drop.intercept <- rep(FALSE, d)
  if (length(drop.intercept) != d) {
    stop("length(drop.intercept) should be equal to number of model formulas")
  }
  lp.overlap <- if (formula$nlp < d) TRUE else FALSE
  G <- gam_setup(
    formula[[1]],
    pterms[[1]],
    data,
    knots,
    sp,
    min.sp,
    H,
    absorb.cons,
    sparse.cons,
    select,
    idLinksBases,
    scale.penalty,
    paraPen,
    gamm.call,
    drop.intercept[1],
    apply.by = apply.by,
    list.call = TRUE,
    modCon = modCon
  )

  G$pterms <- pterms
  G$offset <- list(G$offset)
  G$xlevels <- list(G$xlevels)
  G$assign <- list(G$assign)
  used.sp <- length(G$lsp0)

  if (!is.null(sp) && used.sp > 0) sp <- sp[-(1:used.sp)]
  if (!is.null(min.sp) && nrow(G$L) > 0) min.sp <- min.sp[-(1:nrow(G$L))]

  flpi <- lpi <- list()
  for (i in 1:formula$nlp) lpi[[i]] <- rep(0, 0)
  lpi[[1]] <- 1:ncol(G$X)
  flpi[[1]] <- formula[[1]]$lpi

  pof <- ncol(G$X)
  pstart <- rep(0, d)
  pstart[1] <- 1
  if (d > 1) {
    for (i in 2:d) {
      if (is.null(formula[[i]]$response)) {
        formula[[i]]$response <- formula$response
        mv.response <- FALSE
      } else {
        mv.response <- TRUE
      }
      formula[[i]]$pfok <- 1
      um <- gam_setup(
        formula[[i]],
        pterms[[i]],
        data,
        knots,
        sp,
        min.sp,
        H,
        absorb.cons,
        sparse.cons,
        select,
        idLinksBases,
        scale.penalty,
        paraPen,
        gamm.call,
        drop.intercept[i],
        apply.by = apply.by,
        list.call = TRUE,
        modCon = modCon
      )
      used.sp <- length(um$lsp0)
      if (!is.null(sp) && used.sp > 0) sp <- sp[-(1:used.sp)]
      if (!is.null(min.sp) && nrow(um$L) > 0) min.sp <- min.sp[-(1:nrow(um$L))]

      flpi[[i]] <- formula[[i]]$lpi
      for (j in formula[[i]]$lpi) {
        lpi[[j]] <- c(lpi[[j]], pof + 1:ncol(um$X))
      }
      if (mv.response) G$y <- cbind(G$y, um$y)
      if (i > formula$nlp && !is.null(um$offset)) {
        stop("shared offsets not allowed")
      }
      G$offset[[i]] <- um$offset
      if (!is.null(um$contrasts)) G$contrasts <- c(G$contrasts, um$contrasts)
      G$xlevels[[i]] <- um$xlevels
      G$assign[[i]] <- um$assign
      G$rank <- c(G$rank, um$rank)
      pstart[i] <- pof + 1
      G$X <- cbind(G$X, um$X)
      k <- G$m
      if (um$m) {
        for (j in 1:um$m) {
          um$smooth[[j]]$first.para <- um$smooth[[j]]$first.para + pof
          um$smooth[[j]]$last.para <- um$smooth[[j]]$last.para + pof
          k <- k + 1
          G$smooth[[k]] <- um$smooth[[j]]
        }
      }
      ks <- length(G$S)
      M <- length(um$S)

      if (!is.null(um$L) || !is.null(G$L)) {
        if (is.null(G$L)) G$L <- diag(1, nrow = ks)
        if (is.null(um$L)) um$L <- diag(1, nrow = M)
        G$L <- rbind(
          cbind(G$L, matrix(0, nrow(G$L), ncol(um$L))),
          cbind(matrix(0, nrow(um$L), ncol(G$L)), um$L)
        )
      }

      G$off <- c(G$off, um$off + pof)
      if (M) {
        for (j in 1:M) {
          ks <- ks + 1
          G$S[[ks]] <- um$S[[j]]
        }
      }

      G$m <- G$m + um$m
      G$nsdf[i] <- um$nsdf
      if (!is.null(um$P) || !is.null(G$P)) {
        if (is.null(G$P)) G$P <- diag(1, nrow = pof)
        k <- ncol(um$X)
        if (is.null(um$P)) um$P <- diag(1, nrow = k)
        G$P <- rbind(
          cbind(G$P, matrix(0, pof, k)),
          cbind(matrix(0, k, pof), um$P)
        )
      }
      G$cmX <- c(G$cmX, um$cmX)
      if (um$nsdf > 0) {
        um$term.names[1:um$nsdf] <- paste(
          um$term.names[1:um$nsdf],
          i - 1,
          sep = "."
        )
      }
      G$term.names <- c(G$term.names, um$term.names)
      G$lsp0 <- c(G$lsp0, um$lsp0)
      G$sp <- c(G$sp, um$sp)
      pof <- ncol(G$X)
    }
  }

  if (lp.overlap) {
    rt <- olid(G$X, G$nsdf, pstart, flpi, lpi)
    if (length(rt$dind) > 0) {
      warning(
        "dropping unidentifiable parametric terms from model",
        call. = FALSE
      )
      G$X <- G$X[, -rt$dind]
      G$cmX <- G$cmX[-rt$dind]
      G$term.names <- G$term.names[-rt$dind]
      for (i in 1:length(G$smooth)) {
        k <- sum(rt$dind < G$smooth[[i]]$first.para)
        G$smooth[[i]]$first.para <- G$smooth[[i]]$first.para - k
        G$smooth[[i]]$last.para <- G$smooth[[i]]$last.para - k
      }
      for (i in 1:length(G$off)) G$off[i] <- G$off[i] - sum(rt$dind < G$off[i])
      attr(G$nsdf, "drop.ind") <- rt$dind
    }
  }
  attr(lpi, "overlap") <- lp.overlap
  attr(G$X, "lpi") <- lpi
  attr(G$nsdf, "pstart") <- pstart

  G$g.index <- rep(FALSE, ncol(G$X))
  n.sp0 <- 0
  if (length(G$smooth)) {
    for (i in 1:length(G$smooth)) {
      if (!is.null(G$smooth[[i]]$g.index)) {
        G$g.index[
          G$smooth[[i]]$first.para:G$smooth[[i]]$last.para
        ] <- G$smooth[[i]]$g.index
      }
      n.sp <- length(G$smooth[[i]]$S)
      if (n.sp) {
        G$smooth[[i]]$first.sp <- n.sp0 + 1
        n.sp0 <- G$smooth[[i]]$last.sp <- n.sp0 + n.sp
      }
    }
  }
  if (!any(G$g.index)) G$g.index <- NULL

  G
}


#' Take a set of non-negative integers and return minimal code for generating it
#'
#' This function is derived from \code{mgcv:::compress.iseq}
#'
#' @author Simon N Wood with modifications by Nicholas Clark
#' @noRd
compress_iseq <- function(x) {
  x1 <- sort(x)
  br <- diff(x1) != 1
  txt <- paste(x1[c(TRUE, br)], x1[c(br, TRUE)], sep = ":")
  txt1 <- paste(x1[c(TRUE, br)])
  ii <- x1[c(TRUE, br)] == x1[c(br, TRUE)]
  txt[ii] <- txt1[ii]
  paste("c(", paste(txt, collapse = ","), ")", sep = "")
}

#' Returns a vector dind of columns of X to drop for identifiability
#'
#' This function is derived from \code{mgcv:::olid}
#'
#' @author Simon N Wood with modifications by Nicholas Clark
#' @noRd
olid <- function(X, nsdf, pstart, flpi, lpi) {
  nlp <- length(lpi)
  n <- nrow(X)
  nf <- length(nsdf)
  Xp <- matrix(0, n * nlp, sum(nsdf))
  start <- 1
  ii <- 1:n
  tind <- rep(0, 0)
  for (i in 1:nf) {
    stop <- start - 1 + nsdf[i]
    if (stop >= start) {
      ind <- pstart[i] + 1:nsdf[i] - 1
      for (k in flpi[[i]]) {
        Xp[ii + (k - 1) * n, start:stop] <- X[, ind]
      }
      tind <- c(tind, ind)
      start <- start + nsdf[i]
    }
  }
  qrx <- qr(Xp, LAPACK = TRUE, tol = 0.0)
  r <- mgcv::Rrank(qr.R(qrx))
  if (r == ncol(Xp)) {
    dind <- rep(0, 0)
  } else {
    dind <- tind[sort(qrx$pivot[(r + 1):ncol(X)], decreasing = TRUE)]
    for (d in dind) {
      for (i in 1:nlp) {
        k <- which(d == lpi[[i]])
        if (length(k) > 0) lpi[[i]] <- lpi[[i]][-k]
        k <- which(lpi[[i]] > d)
        if (length(k) > 0) lpi[[i]][k] <- lpi[[i]][k] - 1
      }
    }
  }
  list(dind = dind, lpi = lpi)
}
