#' Define Parameter Bounds for PK Models
#'
#' Utility function to generate lower and upper bounds for
#' pharmacokinetic model parameters, including fixed effects (theta),
#' random effects variances (omega), residual error (sigma),
#' and correlation constraints.
#'
#' @param theta A list with optional elements:
#'   \describe{
#'     \item{lower}{Named list of lower bounds for fixed effects. Defaults to
#'     -Inf for all parameters.}
#'     \item{upper}{Named list of upper bounds for fixed effects. Defaults to
#'     10^9 for all parameters.}
#'   }
#' @param omega A list with optional elements:
#'   \describe{
#'     \item{lower}{Named list of lower bounds for variance terms. Defaults to
#'     10 for all parameters.}
#'     \item{upper}{Named list of upper bounds for variance terms. Defaults to
#'     Inf for all parameters.}
#'   }
#' @param sigma A list with two elements (each itself a list of bounds):
#'   \describe{
#'     \item{add}{Lower and upper bounds for additive error component.
#'     Defaults to 0.001 and Inf.}
#'     \item{prop}{Lower and upper bounds for proportional error component.
#'     Defaults to 0.001 and Inf.}
#'   }
#' @param correlation A list with elements lower and upper
#' giving the bounds for correlation terms. Defaults to 0.1 and 0.8.
#'
#' @returns A named list with four components:
#' \describe{
#'   \item{theta}{List of parameter-specific lower and upper bounds
#'   for fixed effects.}
#'   \item{omega}{List of lower and upper bounds for variance terms.}
#'   \item{sigma}{List with additive (add) and proportional (prop)
#'   error bounds.}
#'   \item{correlation}{List with lower and upper correlation bounds.}
#' }
#'
#' @details
#' Default theta bounds use -Inf for lower limits and 10^9 for upper limits
#' to avoid allowing unrealistically large fixed effect estimates while
#' still providing flexibility during model estimation.
#'
#' @author Zhonghui Huang
#'
#' @examples
#' # Use all default bounds
#'  param.bounds()
#'
#' # Customize only omega lower bounds
#' param.bounds(omega = list(lower = list(cl = 5, vc = 2)))
#'
#' # Adjust sigma proportional error bounds
#' param.bounds(
#'   sigma = list(
#'     add = list(lower = 0.001, upper = 1),
#'     prop = list(lower = 0.01, upper = 0.05)
#'   )
#' )
#'
#' @export
param.bounds <- function(theta = list(lower = NULL, upper = NULL),
                         omega = list(lower = NULL, upper = NULL),
                         sigma = list(add = list(lower = 0.001, upper = Inf),
                                      prop = list(lower = 0.001, upper = Inf)),
                         correlation = list(lower = 0.1, upper = 0.8)) {
  # --- Theta bounds ---
  default.theta.lower <- list(
    ka = -Inf,
    vc = -Inf,
    cl = -Inf,
    vp = -Inf,
    vp2 = -Inf,
    q = -Inf,
    q2 = -Inf,
    tlag = -Inf,
    vmax = -Inf,
    km = -Inf
  )

  default.theta.upper <- list(
    ka = 10^9,
    vc = 10^9,
    cl = 10^9,
    vp = 10^9,
    vp2 = 10^9,
    q = 10^9,
    q2 = 10^9,
    tlag = 10^9,
    vmax = 10^9,
    km = 10^9
  )

  theta.lower <-
    utils::modifyList(default.theta.lower, if (is.null(theta$lower))
      list()
      else
        theta$lower)
  theta.upper <-
    utils::modifyList(default.theta.upper, if (is.null(theta$upper))
      list()
      else
        theta$upper)

  # --- Omega bounds ---
  default.omega.lower <- list(
    ka = 10,
    vc = 10,
    cl = 10,
    vp = 10,
    vp2 = 10,
    q = 10,
    q2 = 10,
    tlag = 10,
    vmax = 10,
    km = 10
  )

  default.omega.upper <- list(
    ka = Inf,
    vc = Inf,
    cl = Inf,
    vp = Inf,
    vp2 = Inf,
    q = Inf,
    q2 = Inf,
    tlag = Inf,
    vmax = Inf,
    km = Inf
  )

  omega.lower <-
    utils::modifyList(default.omega.lower, if (is.null(omega$lower))
      list()
      else
        omega$lower)
  omega.upper <-
    utils::modifyList(default.omega.upper, if (is.null(omega$upper))
      list()
      else
        omega$upper)

  # --- Return ---
  list(
    theta = list(lower = theta.lower, upper = theta.upper),
    omega = list(lower = omega.lower, upper = omega.upper),
    sigma = list(add = sigma$add, prop = sigma$prop),
    correlation = list(lower = correlation$lower, upper = correlation$upper)
  )
}


#' Configure penalty settings for model evaluation
#'
#' Defines rules governing penalty assignment during model adequacy evaluation.
#'
#' @param penalty.value Numeric. Constant penalty assigned to binary violations
#'   and bound constraints.
#'
#' @param step.penalties A named list defining penalty magnitudes used in
#'   step-wise procedures. Each element must contain a numeric vector of length
#'   two representing penalty levels for moderate and critical deviations.
#'
#' @param bounds A list specifying lower and upper parameter limits, as returned
#'   by param.bounds(). The structure can include limits for theta, omega, sigma,
#'   and correlation terms.
#'
#' @param thresholds A named list describing evaluation rules for RSE and
#'   shrinkage. Each component must include a field named method, with value
#'   binary or step, together with the corresponding limit definition:
#'   - If method = binary: a single cutoff value stored in threshold
#'   - If method = step: two deviation boundaries stored in step.levels
#'
#' @param penalty.terms Character vector specifying which components are
#'   considered when penalties are reported. Recognized entries include:
#'   rse, shrinkage, theta, omega, sigma, correlation, covariance, and total.
#'   If total is included, penalties are aggregated across all components and
#'   any other entries are ignored.
#'
#' @details
#' Penalization may be triggered by exceeding predefined parameter bounds
#' (fixed-effect and variance-covariance elements) or by surpassing thresholds
#' for relative standard error (RSE) or shrinkage criteria. Binary and step-wise
#' penalty procedures are supported.
#'
#' @returns A list containing the full penalty configuration for use in fitness().
#'
#' @author Zhonghui Huang
#'
#' @examples
#' # Default configuration
#' penaltyControl()
#'
#' # Custom bounds for selected fixed-effect parameters
#' penaltyControl(bounds = param.bounds(
#'   theta = list(lower = list(cl = 0.01, vc = 0.01))
#' ))
#'
#' # Binary penalty method for RSE
#' penaltyControl(thresholds = list(
#'   rse = list(method = "binary", threshold = 40)
#' ))
#'
#' @seealso \code{\link{param.bounds}()}, \code{\link{fitness}()}.
#'
#' @export

penaltyControl <- function(penalty.value = 10000,
                           step.penalties  = list(
                             rse         = c(10, 10000),
                             shrinkage   = c(10, 10000),
                             bsv         = c(10, 10000),
                             sigma       = list(add  = c(10, 10000),
                                                prop = c(10, 10000)),
                             correlation = c(10, 10000)
                           ),
                           bounds = param.bounds(),
                           thresholds = list(),
                           penalty.terms = c("total")) {
  # Valid penalty terms
  penalty.terms <- tolower(penalty.terms)
  valid.elements <-
    c("rse",
      "shrinkage",
      "theta",
      "omega",
      "sigma",
      "correlation",
      "covariance",
      "total")

  invalid <- setdiff(penalty.terms, valid.elements)
  if (length(invalid) > 0) {
    stop(paste("Invalid penalty.terms:", paste(invalid, collapse = ", ")))
  }
  if ("total" %in% penalty.terms && length(penalty.terms) > 1) {
    warning("'total' overrides other penalty elements - using only 'total'")
    penalty.terms <- "total"
  }

  # Default thresholds for RSE and shrinkage
  default.thresholds <- list(
    rse = list(
      method = "step",
      threshold = 20,
      step.levels = c(20, 40)
    ),
    shrinkage = list(
      method = "step",
      threshold = 30,
      step.levels = c(30, 50)
    ),
    bsv = list(method = "step",
               step.levels = c(10, 10)),
    sigma = list(
      add = list(method = "step",
                 # step.levels = c(0.1, 0.01)
                 step.levels = c(0.001, 0)),
      prop = list(method = "step",
                  # step.levels = c(0.05, 0.01)
                  step.levels = c(0.001, 0))
    ),
    correlation = list(
      method = "step",
      step.levels = list(lower = c(0.1, 0.05),
                         upper = c(0.8, 0.95))
    )
  )

  # Merge user thresholds with defaults
  thresholds <- utils::modifyList(default.thresholds, thresholds)

  # Return full control structure
  return(
    list(
      penalty.value = penalty.value,
      step.penalties = step.penalties,
      bounds = bounds,
      thresholds = thresholds,
      penalty.terms = penalty.terms
    )
  )
}

#' Evaluate fitness of a population pharmacokinetic model
#'
#' Evaluates the quality of a fitted model based on parameter
#' bounds and diagnostic thresholds.
#'
#' @param search.space Character, one of "ivbase" or "oralbase". Default is "ivbase".
#' @param fit Data frame. Model summary from tools such as
#'   `get.mod.lst()`, with parameter estimates and diagnostics.
#' @param dat A data frame containing pharmacokinetic data in standard
#'   nlmixr2 format, including "ID", "TIME", "EVID", and "DV", and may include
#'   additional columns.
#' @param penalty.control List created using \code{penaltyControl()}, including:
#'   \describe{
#'     \item{penalty.value}{Numeric. Default penalty multiplier used in
#'       binary violations.}
#'     \item{step.penalties}{Numeric vector or list. Penalties applied to
#'       step violations (mild, severe).}
#'     \item{bounds}{List of parameter lower/upper bounds, typically from
#'       `param.bounds()`.}
#'     \item{thresholds}{Named list of diagnostic constraints (e.g., RSE,
#'       shrinkage). Each contains a method ("binary" or "step") and the
#'       corresponding threshold or step levels.}
#'     \item{penalty.terms}{Character vector of constraint categories to
#'       penalize. Valid terms include "theta", "rse", "omega",
#'       "shrinkage", "sigma", "correlation", "covariance", and "total".}
#'   }
#' @param objf Character. Column name in fit used as the base objective
#'   function (e.g., "AIC", "BIC", "OBJFV").
#'
#' @returns A data frame extending fit with the following:
#' \itemize{
#'   \item flag.* columns: indicators of constraint violations
#'     (0 = no violation, 1 = mild, 2 = severe).
#'   \item count.constraint.* columns: number of violations per
#'     constraint type.
#'   \item fitness: penalized objective function value, computed from the
#'     specified objf plus applicable penalties.
#' }
#'
#' @author Zhonghui Huang
#'
#' @examples
#' \donttest{
#' # Fit a model (using nlmixr2)
#' pheno <- function() {
#'   ini({
#'     tcl <- log(0.008)
#'     tv <-  log(0.6)
#'     eta.cl + eta.v ~ c(1,
#'                        0.01, 1)
#'     add.err <- 0.1
#'   })
#'   model({
#'     cl <- exp(tcl + eta.cl)
#'     v <- exp(tv + eta.v)
#'     ke <- cl / v
#'     d/dt(A1) = - ke * A1
#'     cp = A1 / v
#'     cp ~ add(add.err)
#'   })
#' }
#' fit <- nlmixr2est::nlmixr2(pheno, pheno_sd, "saem", control = list(print = 0),
#'               table = list(cwres = TRUE, npde = TRUE))
#' Store. <- get.mod.lst(fit.s = fit, 1)
#'  fitness(fit = Store.,dat = pheno_sd)
#' }
#'
#' @seealso \code{\link{penaltyControl}()}, \code{\link{param.bounds}()}.
#' @export

fitness <- function(search.space = "ivbase",
                    fit = NULL,
                    dat = NULL,
                    penalty.control = penaltyControl(),
                    objf = "BIC") {
  # --- Extract penalty controls ---
  penalty.terms     <- penalty.control$penalty.terms
  penalty.value     <- penalty.control$penalty.value
  step.penalties    <- penalty.control$step.penalties
  bounds            <- penalty.control$bounds
  thresholds        <- penalty.control$thresholds

  # --- Thresholds ---
  rse.threshold      <- thresholds$rse$threshold
  rse.method        <- thresholds$rse$method
  rse.steps         <- thresholds$rse$step.levels

  shrinkage.threshold   <- thresholds$shrinkage$threshold
  shrink.method     <- thresholds$shrinkage$method
  shrink.steps <-  thresholds$shrinkage$step.levels

  bsv.method <- thresholds$bsv$method
  bsv.steps  <- thresholds$bsv$step.levels

  sigma.add.method  <- thresholds$sigma$add$method
  sigma.add.steps   <- thresholds$sigma$add$step.levels

  sigma.prop.method <- thresholds$sigma$prop$method
  sigma.prop.steps  <- thresholds$sigma$prop$step.levels

  cor.method <- thresholds$correlation$method
  cor.steps  <- thresholds$correlation$step.levels

  # --- Bounds ---
  theta.lower <- bounds$theta$lower
  theta.upper <- bounds$theta$upper
  bsv.lower   <- bounds$omega$lower
  bsv.upper   <- bounds$omega$upper
  sigma.bounds         <- bounds$sigma
  cor.lower           <- bounds$correlation$lower
  cor.upper           <- bounds$correlation$upper


  colnames(dat) <- toupper(colnames(dat))
  if (!"EVID" %in% colnames(dat))
    stop("dat must contain 'EVID' column.")
  dat.obs <- dat[dat$EVID == 0, ]
  pop.cmin <- stats::aggregate(DV ~ ID, data = dat.obs, FUN = min)
  dat.cmin <- stats::median(pop.cmin$DV)
  sigma.bounds$add$lower <-
    signif(dat.cmin * sigma.bounds$add$lower, 1)

  pop.cmax <- stats::aggregate(DV ~ ID, data = dat.obs, FUN = max)
  dat.cmax <- stats::median(pop.cmax$DV)

  # --- Variable assignments by model type ---
  if (search.space %in% c("ivbase", "oralbase")) {
    param.cols <-
      c(
        "thetacl",
        "thetavc",
        "thetavp",
        "thetavp2",
        "thetaq",
        "thetaq2",
        "thetavmax",
        "thetakm"
      )
    rse.cols <-
      c("rsecl",
        "rsevc",
        "rsevp",
        "rsevp2",
        "rseq",
        "rseq2",
        "rsevmax",
        "rsekm")
    bsv.cols <-
      c("bsvvc",
        "bsvcl",
        "bsvvp",
        "bsvvp2",
        "bsvq",
        "bsvq2",
        "bsvtlag",
        "bsvvmax",
        "bsvkm")
    shrink.cols <-
      c(
        "shrinkcl",
        "shrinkvc",
        "shrinkvp",
        "shrinkvp2",
        "shrinkq",
        "shrinkq2",
        "shrinkvmax",
        "shrinkkm"
      )

    if (search.space == "oralbase") {
      param.cols   <- c("thetaka", param.cols)
      rse.cols     <- c("rseka", rse.cols)
      bsv.cols     <- c("bsvka", bsv.cols)
      shrink.cols  <- c("shrinkka", shrink.cols)
    }
    corr.cols  <- grep("^cor\\.", colnames(fit), value = TRUE)
    error.cols <- c("add", "prop")
  } else if (search.space == "custom") {
    # Set all possible parameter columns for custom search space
    param.cols <- c(
      "thetaka",
      "thetacl",
      "thetavc",
      "thetavp",
      "thetavp2",
      "thetaq",
      "thetaq2",
      "thetavmax",
      "thetakm",
      "thetatlag",
      "thetan",
      "thetamtt",
      "thetabio",
      "thetaD2",
      "thetaF1",
      "thetaFr"
    )

    rse.cols <- c(
      "rseka",
      "rsecl",
      "rsevc",
      "rsevp",
      "rsevp2",
      "rseq",
      "rseq2",
      "rsevmax",
      "rsekm",
      "rsetlag",
      "rsen",
      "rsemtt",
      "rsebio",
      "rseD2",
      "rseF1",
      "rseFr"
    )

    bsv.cols <- c(
      "bsvka",
      "bsvcl",
      "bsvvc",
      "bsvvp",
      "bsvvp2",
      "bsvq",
      "bsvq2",
      "bsvvmax",
      "bsvkm",
      "bsvtlag",
      "bsvn",
      "bsvmtt",
      "bsvbio",
      "bsvD2",
      "bsvF1",
      "bsvFr"
    )

    shrink.cols <- c(
      "shrinkka",
      "shrinkcl",
      "shrinkvc",
      "shrinkvp",
      "shrinkvp2",
      "shrinkq",
      "shrinkq2",
      "shrinkvmax",
      "shrinkkm",
      "shrinktlag",
      "shrinkn",
      "shrinkmtt",
      "shrinkbio",
      "shrinkD2",
      "shrinkF1",
      "shrinkFr"
    )

    corr.cols  <- grep("^cor\\.", colnames(fit), value = TRUE)
    error.cols <- c("add", "prop")

  } else {
    stop(sprintf("Unknown search.space: '%s'", search.space), call. = FALSE)
  }

  # --- Apply parameter bound flags (theta) ---
  for (col in param.cols) {
    flag.col <- paste0("flag.", col)
    fit[[flag.col]] <- 0
    param.name <- sub("theta", "", col)
    fit[[flag.col]][!is.na(fit[[col]]) &
                      fit[[col]] < theta.lower[[param.name]]] <- 1
    fit[[flag.col]][!is.na(fit[[col]]) &
                      fit[[col]] > theta.upper[[param.name]]] <- 1
    if (param.name == "km" && !is.null(dat)) {
      fit[[flag.col]][!is.na(fit[[col]]) & fit[[col]] > dat.cmax] <- 1
    }
  }

  # --- RSE penalty (binary or step) ---
  for (col in rse.cols) {
    flag.col <- paste0("flag.", col)
    fit[[flag.col]] <- 0
    if (rse.method == "binary") {
      fit[[flag.col]][!is.na(fit[[col]]) &
                        fit[[col]] > rse.threshold] <- 1
    } else if (rse.method == "step" && !is.null(rse.steps)) {
      fit[[flag.col]][!is.na(fit[[col]]) &
                        fit[[col]] > rse.steps[1] &
                        fit[[col]] <= rse.steps[2]] <- 1
      fit[[flag.col]][!is.na(fit[[col]]) &
                        fit[[col]] > rse.steps[2]] <- 2
    }
  }

  # --- Omega (BSV) bounds ---
  for (col in bsv.cols) {
    flag.col <- paste0("flag.", col)
    fit[[flag.col]] <- 0
    value <- fit[[col]]

    if (bsv.method == "binary") {
      fit[[flag.col]][!is.na(value) &
                        value < bsv.lower[[sub("bsv", "", col)]]] <-
        1
    } else if (bsv.method == "step" && !is.null(bsv.steps)) {
      # Severe: value < step1
      fit[[flag.col]][!is.na(value) & value < bsv.steps[2]] <- 2
      # Mild: step1 ≤ value < step2
      fit[[flag.col]][!is.na(value) &
                        value >= bsv.steps[2] &
                        value < bsv.steps[1]] <- 1
    }
  }


  # --- Sigma (residual error) ---
  for (col in error.cols) {
    flag.col <- paste0("flag.", col)
    fit[[flag.col]] <- 0
    value <- fit[[col]]

    # Additive Error (add)
    if (col == "add") {
      if (sigma.add.method == "step" && !is.null(sigma.add.steps)) {
        fit[[flag.col]][!is.na(value) &
                          value <= sigma.add.steps[1] &
                          value > sigma.add.steps[2]] <- 1
        fit[[flag.col]][!is.na(value) &
                          value < sigma.add.steps[2]] <- 2
      } else if (sigma.add.method == "binary") {
        fit[[flag.col]][!is.na(value) & value < sigma.bounds$add$lower] <- 1
      }
    }

    # Proportional Error (prop)
    if (col == "prop") {
      if (sigma.prop.method == "step" && !is.null(sigma.prop.steps)) {
        fit[[flag.col]][!is.na(value) &
                          value <= sigma.prop.steps[1] &
                          value > sigma.prop.steps[2]] <- 1
        fit[[flag.col]][!is.na(value) &
                          value < sigma.prop.steps[2]] <- 2
      } else if (sigma.prop.method == "binary") {
        fit[[flag.col]][!is.na(value) &
                          value < sigma.bounds$prop$lower] <- 1
      }
    }
  }


  # --- Shrinkage (binary or step) ---
  for (col in shrink.cols) {
    flag.col <- paste0("flag.", col)
    fit[[flag.col]] <- 0
    if (shrink.method == "binary") {
      fit[[flag.col]][!is.na(fit[[col]]) &
                        fit[[col]] > shrinkage.threshold] <- 1
    } else if (shrink.method == "step" && !is.null(shrink.steps)) {
      fit[[flag.col]][!is.na(fit[[col]]) &
                        fit[[col]] > shrink.steps[1] &
                        fit[[col]] <= shrink.steps[2]] <- 1
      fit[[flag.col]][!is.na(fit[[col]]) &
                        fit[[col]] > shrink.steps[2]] <- 2
    }
  }

  # --- Correlation ---
  for (col in corr.cols) {
    flag.col <- paste0("flag.", col)
    fit[[flag.col]] <- 0
    value <- fit[[col]]
    if (cor.method == "binary") {
      fit[[flag.col]][!is.na(value) & value != 0 &
                        (abs(value) > cor.upper |
                           abs(value) < cor.lower)] <- 1
    } else if (cor.method == "step" && !is.null(cor.steps)) {
      value.abs <- abs(value)
      # ---- Lower bound ----
      fit[[flag.col]][!is.na(value.abs) &
                        value.abs < cor.steps$lower[1] &
                        value.abs >= cor.steps$lower[2]] <- 1
      fit[[flag.col]][!is.na(value.abs) &
                        value.abs < cor.steps$lower[2]] <- 2
      # ---- Upper bound ----
      fit[[flag.col]][!is.na(value.abs) &
                        value.abs > cor.steps$upper[1] &
                        value.abs < cor.steps$upper[2]] <- 1
      fit[[flag.col]][!is.na(value.abs) &
                        value.abs >= cor.steps$upper[2] &
                        value.abs < 1] <- 2
      # ---- Perfect correlation ----
      fit[[flag.col]][!is.na(value.abs) & value.abs == 1] <- 2
      # ---- Zero correlation ----
      fit[[flag.col]][value.abs == 0 | is.na(value.abs)] <- 0
    }
  }

  # --- Covariance method flag ---
  fit$flag.covariance <- 1
  fit$flag.covariance[fit$model.covMethod %in% c("r", "s", "r,s", "linFim")] <-
    0

  # --- Collect flag columns ---
  cols.theta   <- grep("^flag\\.thet", colnames(fit), value = TRUE)
  cols.rse     <- grep("^flag\\.rse", colnames(fit), value = TRUE)
  cols.omega   <- grep("^flag\\.bsv", colnames(fit), value = TRUE)
  cols.shrink  <-
    grep("^flag\\.shrink", colnames(fit), value = TRUE)
  cols.corr    <- grep("^flag\\.cor", colnames(fit), value = TRUE)
  cols.sigma   <-
    grep("^flag\\.(add|prop)$", colnames(fit), value = TRUE)
  flag.cols    <- grep("^flag\\.", colnames(fit), value = TRUE)

  # --- Count constraint violations ---
  fit$count.constraint.theta <-
    if (length(cols.theta) > 0) {
      rowSums(fit[, cols.theta, drop = FALSE] > 0, na.rm = TRUE)
    } else {
      rep(0, nrow(fit))
    }

  fit$count.constraint.rse <-
    if (length(cols.rse) > 0) {
      rowSums(fit[, cols.rse, drop = FALSE] > 0, na.rm = TRUE)
    } else {
      rep(0, nrow(fit))
    }

  fit$count.constraint.omega <-
    if (length(cols.omega) > 0) {
      rowSums(fit[, cols.omega, drop = FALSE] > 0, na.rm = TRUE)
    } else {
      rep(0, nrow(fit))
    }

  fit$count.constraint.shrinkage <-
    if (length(cols.shrink) > 0) {
      rowSums(fit[, cols.shrink, drop = FALSE] > 0, na.rm = TRUE)
    } else {
      rep(0, nrow(fit))
    }

  fit$count.constraint.correlation <-
    if (length(cols.corr) > 0) {
      rowSums(fit[, cols.corr, drop = FALSE] > 0, na.rm = TRUE)
    } else {
      rep(0, nrow(fit))
    }

  fit$count.constraint.sigma <-
    if (length(cols.sigma) > 0) {
      rowSums(fit[, cols.sigma, drop = FALSE] > 0, na.rm = TRUE)
    } else {
      rep(0, nrow(fit))
    }

  fit$count.constraint.total <-
    if (length(flag.cols) > 0) {
      rowSums(fit[, flag.cols, drop = FALSE] > 0, na.rm = TRUE)
    } else {
      rep(0, nrow(fit))
    }

  # --- Fitness calculation ---
  valid.objectives <- c("AIC", "BIC", "OBJFV")
  if (!(objf %in% valid.objectives))
    stop(sprintf("Invalid objective '%s'.", objf))
  if (!(objf %in% colnames(fit)))
    stop(sprintf("Column '%s' not found in fit.", objf))

  fit$fitness <- fit[[objf]]

  if ("total" %in% penalty.terms) {
    penalty.terms <- c("theta",
                       "omega",
                       "correlation",
                       "sigma",
                       "covariance",
                       "rse",
                       "shrinkage")
  }
  if ("theta" %in% penalty.terms) {
    fit$fitness <-
      fit$fitness + penalty.value * fit$count.constraint.theta
  }
  # if ("omega" %in% penalty.terms) {
  #   fit$fitness <-
  #     fit$fitness + penalty.value * fit$count.constraint.omega
  # }
  if ("covariance" %in% penalty.terms) {
    fit$fitness <- fit$fitness + penalty.value * fit$flag.covariance
  }

  # Initialize a penalty accumulator (one value per row in the dataset)
  penalty_acc <- numeric(nrow(fit))
  # ---- RSE Penalty ----
  if ("rse" %in% penalty.terms) {
    for (col in cols.rse) {
      if (rse.method == "binary") {
        # Add penalty if RSE value is 1 or NA (indicating estimation failure or boundary)
        penalty_acc <- penalty_acc +
          penalty.value * ifelse(is.na(fit[[col]]) |
                                   fit[[col]] == 1, 1, 0)
      } else if (rse.method == "step") {
        penalty_acc <- penalty_acc +
          step.penalties$rse[1] * ifelse(is.na(fit[[col]]) |
                                           fit[[col]] == 1, 1, 0) +
          step.penalties$rse[2] * ifelse(fit[[col]] == 2, 1, 0)
      }
    }
  }

  # ---- SHRINKAGE Penalty ----
  if ("shrinkage" %in% penalty.terms) {
    for (col in cols.shrink) {
      if (shrink.method == "binary") {
        # Add penalty if shrinkage flag is 1 or NA
        penalty_acc <- penalty_acc +
          penalty.value * ifelse(is.na(fit[[col]]) |
                                   fit[[col]] == 1, 1, 0)
      } else if (shrink.method == "step") {
        penalty_acc <- penalty_acc +
          step.penalties$shrinkage[1] * ifelse(is.na(fit[[col]]) |
                                                 fit[[col]] == 1, 1, 0) +
          step.penalties$shrinkage[2] * ifelse(fit[[col]] == 2, 1, 0)
      }
    }
  }

  # Sigma Penalties for add and prop residual errors
  if ("sigma" %in% penalty.terms) {
    for (col in cols.sigma) {
      if (col == "flag.add") {
        if (sigma.add.method == "binary") {
          penalty_acc <- penalty_acc +
            penalty.value * ifelse(is.na(fit[[col]]) |
                                     fit[[col]] == 1, 1, 0)
        } else if (sigma.add.method == "step") {
          penalty_acc <- penalty_acc +
            step.penalties$sigma$add[1] * ifelse(is.na(fit[[col]]) |
                                                   fit[[col]] == 1, 1, 0) +
            step.penalties$sigma$add[2] * ifelse(fit[[col]] == 2, 1, 0)
        }
      }
      if (col == "flag.prop") {
        if (sigma.prop.method == "binary") {
          penalty_acc <- penalty_acc +
            penalty.value * ifelse(is.na(fit[[col]]) |
                                     fit[[col]] == 1, 1, 0)
        } else if (sigma.prop.method == "step") {
          penalty_acc <- penalty_acc +
            step.penalties$sigma$prop[1] * ifelse(is.na(fit[[col]]) |
                                                    fit[[col]] == 1, 1, 0) +
            step.penalties$sigma$prop[2] * ifelse(fit[[col]] == 2, 1, 0)
        }
      }
    }
  }

  if ("omega" %in% penalty.terms) {
    for (col in cols.omega) {
      if (bsv.method == "binary") {
        penalty_acc <- penalty_acc +
          penalty.value * ifelse(is.na(fit[[col]]) |
                                   fit[[col]] == 1, 1, 0)
      } else if (bsv.method == "step") {
        penalty_acc <- penalty_acc +
          step.penalties$bsv[1] * ifelse(is.na(fit[[col]]) |
                                           fit[[col]] == 1, 1, 0) +
          step.penalties$bsv[2] * ifelse(fit[[col]] == 2, 1, 0)
      }
    }
  }

  if ("correlation" %in% penalty.terms) {
    for (col in cols.corr) {
      if (cor.method == "binary") {
        penalty_acc <-
          penalty_acc + penalty.value * ifelse(is.na(fit[[col]]) |
                                                 fit[[col]] == 1, 1, 0)
      } else if (cor.method == "step") {
        penalty_acc <- penalty_acc +
          step.penalties$correlation[1] * ifelse(is.na(fit[[col]]) |
                                                   fit[[col]] == 1, 1, 0) +
          step.penalties$correlation[2] * ifelse(fit[[col]] == 2, 1, 0)
      }
    }
  }
  # ---- Final fitness update ----
  # Add the total accumulated penalty to the fitness for each model
  fit$fitness <- fit$fitness + penalty_acc

  return(fit)
}
