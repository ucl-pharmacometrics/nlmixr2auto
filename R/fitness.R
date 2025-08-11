#' Construct default parameter bounds for model fitness evaluation
#'
#' @param theta List of lower and upper bounds for fixed-effect parameters (theta)
#' @param omega List of lower bounds for BSV parameters (omega)
#' @param sigma List with named lower bounds for residual error (additive and proportional)
#' @param correlation List of lower and upper bounds for interindividual correlation
#'
#' @return A structured list containing parameter bounds
#' @export
param.bounds <- function(theta = list(lower = NULL, upper = NULL),
                         omega = list(lower = NULL, upper = NULL),
                           sigma = list(add = list(lower = NULL, upper = Inf),
                                      prop = list(lower = 0.05, upper = Inf)),
                         correlation = list(lower = 0.1, upper = Inf)) {
  # --- Theta bounds ---
  default.theta.lower <- list(
    ka = 0.001,
    vc = 0.001,
    cl = 0.001,
    vp = 0.001,
    vp2 = 0.001,
    q = 0.001,
    q2 = 0.001,
    tlag = 0.001,
    vmax = 0.001,
    km = 0.001
  )
  
  default.theta.upper <- list(
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
  
  theta.lower <-
    modifyList(default.theta.lower, if (is.null(theta$lower))
      list()
      else
        theta$lower)
  theta.upper <-
    modifyList(default.theta.upper, if (is.null(theta$upper))
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
    modifyList(default.omega.lower, if (is.null(omega$lower))
      list()
      else
        omega$lower)
  omega.upper <-
    modifyList(default.omega.upper, if (is.null(omega$upper))
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


#' Define penalty control settings for model fitness evaluation
#'
#' This function sets up penalty rules used in evaluating model quality based on
#' bounds (e.g., parameter lower/upper limits) and threshold-based metrics such as
#' RSE and shrinkage. Supports both binary and stepwise penalty schemes.
#'
#' @param penalty.value Numeric. Default penalty applied for binary violations and bound constraints.
#' @param step.penalties Numeric vector of length 2. Penalties for mild and severe violations
#'        in step-based methods (e.g., c(mild = 10, severe = 10000)).
#' @param bounds A list of parameter bounds defined by `param.bounds()` (includes theta, omega, sigma, correlation).
#' @param thresholds A named list of threshold control for RSE and shrinkage.
#'        Each element should contain a `method` ("binary" or "step"), and:
#'        - For "binary": `threshold` (single value)
#'        - For "step": `step.levels` (e.g., c(30, 50))
#' @param penalty.terms Character vector of constraint types to penalize.
#'        Valid terms: "rse", "shrinkage", "theta", "omega", "sigma",
#'        "correlation", "covariance", "total".
#'
#' @return A list containing penalty configuration for use in the `fitness()` function.
#'
#' @examples
#' # Use all defaults
#' penaltyControl()
#'
#' # Override only theta bounds
#' penaltyControl(bounds = param.bounds(
#'   theta = list(lower = list(cl = 0.01, vc = 0.01))
#' ))
#'
#' # Use binary method for RSE with a custom threshold
#' penaltyControl(thresholds = list(
#'   rse = list(method = "binary", threshold = 40)
#' ))
#'
#'
#' @export
penaltyControl <- function(penalty.value = 10000,
                           step.penalties = c(10, 10000),
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
    warning("'total' overrides other penalty elements — using only 'total'")
    penalty.terms <- "total"
  }
  
  # Default thresholds for RSE and shrinkage
  default.thresholds <- list(
    rse = list(
      method = "step",
      threshold = 20,
      step.levels = c(20, 30)
    ),
    shrinkage = list(
      method = "step",
      threshold = 30,
      step.levels = c(30, 50)
    ),
    bsv = list(
      method = "step",
      step.levels = c(10, 5)  # step1 = severe threshold, step2 = mild threshold
    ),
    sigma = list(
      add = list(
        method = "step",       
        step.levels = c(0.01, 0) 
      ),
      prop = list(
        method = "step",       
        step.levels = c(0.05, 0) 
      )
    ),
    correlation = list(
      method = "step",
      step.levels = list(lower = 0.1, upper = 0.8)
    )
  )
  
  # Merge user thresholds with defaults
  thresholds <- modifyList(default.thresholds, thresholds)
  
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
#' This function evaluates the quality of a fitted model based on parameter bounds and
#' diagnostic thresholds (e.g., relative standard error [RSE], shrinkage, residual error,
#' and inter-individual correlation). Violations are flagged and penalized, either via
#' binary (yes/no) or stepwise (graded) penalties.
#'
#' @param search.space Character. Defines the structural model type (e.g., `"ivbase"`, `"oralbase"`).
#' @param fit Data frame. Model summary from tools such as `get.mod.lst()`, with parameter estimates and diagnostics.
#' @param dat Data frame. Observed data containing columns `EVID`, `DV`, and `ID`; used for estimating default additive error if missing.
#' @param penalty.control List. Created using `penaltyControl()`, includes:
#'   \describe{
#'     \item{\code{penalty.value}}{Numeric. Default penalty multiplier used in binary violations.}
#'     \item{\code{step.penalties}}{Numeric vector of length 2. Penalties applied to step violations (mild, severe).}
#'     \item{\code{bounds}}{List. Parameter lower/upper bounds, created via `param.bounds()`.}
#'     \item{\code{thresholds}}{Named list. Diagnostic constraints (e.g., RSE, shrinkage), each with `method` ("binary" or "step"), and corresponding threshold or step levels.}
#'     \item{\code{penalty.terms}}{Character vector. Constraint categories to penalize. Valid terms include `"theta"`, `"rse"`, `"omega"`, `"shrinkage"`, `"sigma"`, `"correlation"`, `"covariance"`, and `"total"`.}
#'   }
#' @param objf Character. Column name in `fit` used as the base objective function (e.g., `"AIC"`, `"BIC"`, `"OBJFV"`).
#'
#' @return A data frame extending `fit` with the following:
#' \itemize{
#'   \item \code{flag.*} columns: indicators of constraint violations (0 = no violation, 1 = mild, 2 = severe).
#'   \item \code{count.constraint.*} columns: number of violations per constraint type.
#'   \item \code{fitness}: penalized objective function value, computed from the specified `objf` plus applicable penalties.
#' }
#'
#' @examples
#' \dontrun{
#' # Define the model
#' pheno <- function() {
#'   ini({
#'     lcl <- log(0.008) # typical value of clearance
#'     lvc  <- log(0.6)   # typical value of volume
#'     eta.cl + eta.vc ~ c(1,
#'                        0.01, 1)
#'     add.err <- 0.1    # residual variability
#'   })
#'   model({
#'     cl <- exp(lcl + eta.cl) # individual value of clearance
#'     vc  <- exp(lvc + eta.vc)   # individual value of volume
#'     ke <- cl / vc            # elimination rate constant
#'     d/dt(A1) = - ke * A1    # model differential equation
#'     cp = A1 / vc             # concentration in plasma
#'     cp ~ add(add.err)       # define error model
#'   })
#' }

#' fit <- nlmixr(pheno, pheno_sd, "saem", control = list(print = 0),
#'               table = list(cwres = TRUE, npde = TRUE))

#' Store. <- get.mod.lst(fit.s = fit, 1)
#' # --- Example 1: Default settings (step penalties) ---
#' result.default <- fitness(
#'   fit = Store.,
#'   dat = pheno_sd
#' )
#' print(result.default)
#' # --- Example 2: Use binary penalties for RSE and shrinkage ---
#' result.binary <- fitness(
#'   fit = Store.,
#'   dat = pheno_sd,
#'   penalty.control = penaltyControl(
#'     thresholds = list(
#'       rse = list(method = "binary"),
#'       shrinkage = list(method = "binary")
#'     )
#'   )
#' )
#' print(result.binary)
#' }
#' 
#' @export
fitness <- function(search.space = "ivbase",
                    fit = NULL,
                    dat = NULL,
                    penalty.control = penaltyControl(),
                    objf = "AIC") {
  
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
  
  # --- Precompute cadd if missing ---
  if (!is.null(dat)) {
    colnames(dat) <- toupper(colnames(dat))
    if (!"EVID" %in% colnames(dat)) stop("dat must contain 'EVID' column.")
    dat.obs <- dat[dat$EVID == 0, ]
    pop.cmax <- aggregate(DV ~ ID, data = dat.obs, FUN = max)
    dat.cmax <- median(pop.cmax$DV)
  }
  
  if (is.null(sigma.bounds$add$lower)) {
    sigma.bounds$add$lower <- round(dat.cmax / 10000, 0)
  }
  
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
  }
  
  # --- Apply parameter bound flags (theta) ---
  for (col in param.cols) {
    flag.col <- paste0("flag.", col)
    fit[[flag.col]] <- 0
    param.name <- sub("theta", "", col)
    fit[[flag.col]][!is.na(fit[[col]]) & fit[[col]] < theta.lower[[param.name]]] <- 1
    fit[[flag.col]][!is.na(fit[[col]]) & fit[[col]] > theta.upper[[param.name]]] <- 1
    if (param.name == "km" && !is.null(dat)) {
      fit[[flag.col]][!is.na(fit[[col]]) & fit[[col]] > dat.cmax] <- 1
    }
  }
  
  # --- RSE penalty (binary or step) ---
  for (col in rse.cols) {
    flag.col <- paste0("flag.", col)
    fit[[flag.col]] <- 0
    if (rse.method == "binary") {
      fit[[flag.col]][!is.na(fit[[col]]) & fit[[col]] > rse.threshold] <- 1
    } else if (rse.method == "step" && !is.null(rse.steps)) {
      fit[[flag.col]][!is.na(fit[[col]]) & fit[[col]] > rse.steps[1] & fit[[col]] <= rse.steps[2]] <- 1
      fit[[flag.col]][!is.na(fit[[col]]) & fit[[col]] > rse.steps[2]] <- 2
    }
  }
  
  # --- Omega (BSV) bounds ---
  for (col in bsv.cols) {
    flag.col <- paste0("flag.", col)
    fit[[flag.col]] <- 0
    value <- fit[[col]]
    
    if (bsv.method == "binary") {
      fit[[flag.col]][!is.na(value) & value < bsv.lower[[sub("bsv", "", col)]]] <- 1
    } else if (bsv.method == "step" && !is.null(bsv.steps)) {
      # Severe: value < step1
      fit[[flag.col]][!is.na(value) & value < bsv.steps[2]] <- 2
      # Mild: step1 ≤ value < step2
      fit[[flag.col]][!is.na(value) & value >= bsv.steps[2] & value < bsv.steps[1]] <- 1
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
        fit[[flag.col]][!is.na(value) & value <= sigma.add.steps[1] & value > sigma.add.steps[2]] <- 1
        fit[[flag.col]][!is.na(value) & value < sigma.add.steps[2]] <- 2
      } else if (sigma.add.method == "binary" ) {
        fit[[flag.col]][!is.na(value) & value < sigma.bounds$add$lower] <- 1
      }
    }
    
    # Proportional Error (prop)
    if (col == "prop") {
      if (sigma.prop.method == "step" && !is.null(sigma.prop.steps)) {
        fit[[flag.col]][!is.na(value) & value <= sigma.prop.steps[1] & value > sigma.prop.steps[2]] <- 1
        fit[[flag.col]][!is.na(value) & value < sigma.prop.steps[2]] <- 2
      } else if (sigma.prop.method == "binary") {
        fit[[flag.col]][!is.na(value) & value < sigma.bounds$prop$lower] <- 1
      }
    }
  }
  
  
  # --- Shrinkage (binary or step) ---
  for (col in shrink.cols) {
    flag.col <- paste0("flag.", col)
    fit[[flag.col]] <- 0
    if (shrink.method == "binary") {
      fit[[flag.col]][!is.na(fit[[col]]) & fit[[col]] > shrinkage.threshold] <- 1
    } else if (shrink.method == "step" && !is.null(shrink.steps)) {
      fit[[flag.col]][!is.na(fit[[col]]) &
                        fit[[col]] > shrink.steps[1] & fit[[col]] <= shrink.steps[2]] <- 1
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
                        (abs(value) >= cor.steps$upper | abs(value) <= cor.steps$lower)] <- 1
    } else if (cor.method == "step" && !is.null(cor.steps)) {
      if (cor.steps$lower >= cor.steps$upper) stop("Correlation step thresholds invalid: lower >= upper")
      fit[[flag.col]][!is.na(value) & abs(value) == 1] <- 2
      fit[[flag.col]][!is.na(value) & abs(value) > 0 & abs(value) <= cor.steps$lower] <- 1
      fit[[flag.col]][!is.na(value) & abs(value) >= cor.steps$upper & abs(value) < 1] <- 1
      fit[[flag.col]][abs(value) == 0 | is.na(value)] <- 0
    }
    
  }
  
  
  # --- Covariance method flag ---
  fit$flag.covariance <- 1
  fit$flag.covariance[fit$model.covMethod %in% c("r", "s", "r,s", "linFim")] <- 0
  
  # --- Collect flag columns ---
  cols.theta   <- grep("^flag\\.thet", colnames(fit), value = TRUE)
  cols.rse     <- grep("^flag\\.rse", colnames(fit), value = TRUE)
  cols.omega   <- grep("^flag\\.bsv", colnames(fit), value = TRUE)
  cols.shrink  <- grep("^flag\\.shrink", colnames(fit), value = TRUE)
  cols.corr    <- grep("^flag\\.cor", colnames(fit), value = TRUE)
  cols.sigma   <- grep("^flag\\.(add|prop)$", colnames(fit), value = TRUE)
  flag.cols    <- grep("^flag\\.", colnames(fit), value = TRUE)
  
  # --- Count constraint violations ---
  fit$count.constraint.theta <-
    if (length(cols.theta) > 0) {
      rowSums(fit[, cols.theta, drop = FALSE])
    } else {
      rep(0, nrow(fit))
    }
  
  fit$count.constraint.rse <-
    if (length(cols.rse) > 0) {
      rowSums(fit[, cols.rse, drop = FALSE])
    } else {
      rep(0, nrow(fit))
    }
  
  fit$count.constraint.omega <-
    if (length(cols.omega) > 0) {
      rowSums(fit[, cols.omega, drop = FALSE])
    } else {
      rep(0, nrow(fit))
    }
  
  fit$count.constraint.shrinkage <-
    if (length(cols.shrink) > 0) {
      rowSums(fit[, cols.shrink, drop = FALSE])
    } else {
      rep(0, nrow(fit))
    }
  
  fit$count.constraint.correlation <-
    if (length(cols.corr) > 0) {
      rowSums(fit[, cols.corr, drop = FALSE])
    } else {
      rep(0, nrow(fit))
    }
  
  fit$count.constraint.sigma <-
    if (length(cols.sigma) > 0) {
      rowSums(fit[, cols.sigma, drop = FALSE])
    } else {
      rep(0, nrow(fit))
    }
  
  fit$count.constraint.total <-
    if (length(flag.cols) > 0) {
      rowSums(fit[, flag.cols, drop = FALSE])
    } else {
      rep(0, nrow(fit))
    }
  
  # --- Fitness calculation ---
  valid.objectives <- c("AIC", "BIC", "OBJFV")
  if (!(objf %in% valid.objectives)) stop(sprintf("Invalid objective '%s'.", objf))
  if (!(objf %in% colnames(fit))) stop(sprintf("Column '%s' not found in fit.", objf))
  
  fit$fitness <- fit[[objf]]
  
    if ("total" %in% penalty.terms) {
      penalty.terms <- c(
        "theta", "omega", "correlation", 
        "sigma", "covariance", "rse", "shrinkage"
      )
    }
    if ("theta" %in% penalty.terms) {
      fit$fitness <-
        fit$fitness + penalty.value * fit$count.constraint.theta
    }
    if ("omega" %in% penalty.terms) {
      fit$fitness <-
        fit$fitness + penalty.value * fit$count.constraint.omega
    }
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
          penalty.value * ifelse(is.na(fit[[col]]) | fit[[col]] == 1, 1, 0)
      } else if (rse.method == "step") {
        # Apply stepwise penalties for RSE values of 1 and 2
        penalty_acc <- penalty_acc +
          step.penalties[1] * ifelse(is.na(fit[[col]]) | fit[[col]] == 1, 1, 0) +
          step.penalties[2] * ifelse(fit[[col]] == 2, 1, 0)
      }
    }
  }
  
  # ---- SHRINKAGE Penalty ----
  if ("shrinkage" %in% penalty.terms) {
    for (col in cols.shrink) {
      if (shrink.method == "binary") {
        # Add penalty if shrinkage flag is 1 or NA
        penalty_acc <- penalty_acc +
          penalty.value * ifelse(is.na(fit[[col]]) | fit[[col]] == 1, 1, 0)
      } else if (shrink.method == "step") {
        # Apply stepwise penalties for shrinkage values of 1 and 2
        penalty_acc <- penalty_acc +
          step.penalties[1] * ifelse(is.na(fit[[col]]) | fit[[col]] == 1, 1, 0) +
          step.penalties[2] * ifelse(fit[[col]] == 2, 1, 0)
      }
    }
  }
  
  # Sigma Penalties for add and prop residual errors
  if ("sigma" %in% penalty.terms) {
    for (col in cols.sigma) {
      if (col == "flag.add") {
        if (sigma.add.method == "binary") {
          penalty_acc <- penalty_acc +
            penalty.value * ifelse(is.na(fit[[col]]) | fit[[col]] == 1, 1, 0)
        } else if (sigma.add.method == "step") {
          penalty_acc <- penalty_acc +
            step.penalties[1] * ifelse(is.na(fit[[col]]) | fit[[col]] == 1, 1, 0) +
            step.penalties[2] * ifelse(fit[[col]] == 2, 1, 0)
        }
      }
      if (col == "flag.prop") {
        if (sigma.prop.method == "binary") {
          penalty_acc <- penalty_acc +
            penalty.value * ifelse(is.na(fit[[col]]) | fit[[col]] == 1, 1, 0)
        } else if (sigma.prop.method == "step") {
          penalty_acc <- penalty_acc +
            step.penalties[1] * ifelse(is.na(fit[[col]]) | fit[[col]] == 1, 1, 0) +
            step.penalties[2] * ifelse(fit[[col]] == 2, 1, 0)
        }
      }
    }
  }
  
  if ("correlation" %in% penalty.terms) {
    for (col in cols.corr) {
      if (cor.method == "binary") {
        penalty_acc <- penalty_acc + penalty.value * ifelse(is.na(fit[[col]]) | fit[[col]] == 1, 1, 0)
      } else if (cor.method == "step") {
        penalty_acc <- penalty_acc +
          step.penalties[1] * ifelse(is.na(fit[[col]]) | fit[[col]] == 1, 1, 0) +
          step.penalties[2] * ifelse(fit[[col]] == 2, 1, 0)
      }
    }
  }
  # ---- Final fitness update ----
  # Add the total accumulated penalty to the fitness for each model
  fit$fitness <- fit$fitness + penalty_acc

  return(fit)
}

