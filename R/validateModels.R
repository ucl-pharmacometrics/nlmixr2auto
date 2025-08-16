#' Validate and Correct Model Codes
#'
#' @description
#' This function validates and corrects model codes for pharmacokinetic model structures.
#' It is designed to handle inputs from Genetic Algorithm (GA), Ant Colony Optimization (ACO),
#' and Tabu Search (TS) workflows, ensuring that only valid model structures are retained.
#'
#' ## Processing Logic
#' - **GA**:
#'   1. Input: GA binary code.
#'   2. Convert to categorical code via \code{\link{parseCode}}.
#'   3. Apply model validity corrections (compartment logic, Michaelis–Menten, oral absorption).
#'   4. Convert corrected categorical code back to GA binary code
#'      (no.cmpt1/2, rv1/2 encoding).
#'
#' - **ACO** and **TS**:
#'   1. Input: categorical code.
#'   2. Apply model validity corrections directly.
#'   3. Output corrected categorical code.
#'
#' @param string Numeric vector representing either a GA binary code
#'   or a categorical code (ACO/TS).
#' @param search.space Character: `"ivbase"` or `"oralbase"`.
#'   Determines the model structure space (intravenous vs oral).
#' @param code.source Character: `"GA"`, `"ACO"`, or `"TS"`.
#'   Specifies the origin of the code:
#'   - `"GA"` = binary code,
#'   - `"ACO"` or `"TS"` = categorical code.
#'
#' @return A numeric vector:
#'   - For GA input: corrected GA binary code.
#'   - For ACO/TS input: corrected categorical code.
#'
#' @seealso \code{\link{parseCode}}
#' @export
#'
validateModels <-
  function(string,
           search.space = "ivbase",
           code.source = "GA") {
    # --- Step 1: If GA, parse binary code into categorical code ---
    if (code.source == "GA") {
      string <- parseCode(string, search.space)
    }

    # --- Step 2: Extract variables according to search space ---
    if (search.space == "ivbase") {
      no.cmpt <- string[1]
      eta.km  <- string[2]
      eta.vc  <- string[3]
      eta.vp  <- string[4]
      eta.vp2 <- string[5]
      eta.q   <- string[6]
      eta.q2  <- string[7]
      mm      <- string[8]
      mcorr   <- string[9]
      rv      <- string[10]
      eta.ka  <- NULL
    } else if (search.space == "oralbase") {
      no.cmpt <- string[1]
      eta.km  <- string[2]
      eta.vc  <- string[3]
      eta.vp  <- string[4]
      eta.vp2 <- string[5]
      eta.q   <- string[6]
      eta.q2  <- string[7]
      eta.ka  <- string[8]
      mm      <- string[9]
      mcorr   <- string[10]
      rv      <- string[11]
    }

    if (mm==0){
    # Build base data frame
    string.df <- data.frame(
      no.cmpt = no.cmpt,
      eta.vmax=0,
      eta.km = eta.km,
      eta.cl = 1,
      eta.vc = eta.vc,
      eta.vp = eta.vp,
      eta.vp2 = eta.vp2,
      eta.q = eta.q,
      eta.q2 = eta.q2,
      mm = mm,
      mcorr = mcorr,
      rv = rv,
      stringsAsFactors = FALSE
    )
    }

    if (mm==1){
      # Build base data frame
      string.df <- data.frame(
        no.cmpt = no.cmpt,
        eta.vmax=1,
        eta.km = eta.km,
        eta.cl = 0,
        eta.vc = eta.vc,
        eta.vp = eta.vp,
        eta.vp2 = eta.vp2,
        eta.q = eta.q,
        eta.q2 = eta.q2,
        mm = mm,
        mcorr = mcorr,
        rv = rv,
        stringsAsFactors = FALSE
      )
    }

    # If oral, insert eta.ka
    if (!is.null(eta.ka)) {
      pos <- which(names(string.df) == "eta.q2")
      string.df <- cbind(string.df[, 1:pos, drop = FALSE],
                         eta.ka = eta.ka,
                         string.df[, (pos + 1):ncol(string.df), drop = FALSE])
    }

    # --- Step 3: Apply corrections ---
    # Compartment-based corrections
    if (all(c("no.cmpt", "eta.vp") %in% names(string.df)))
      string.df$eta.vp [string.df$no.cmpt == 1] <- 0
    if (all(c("no.cmpt", "eta.q") %in% names(string.df)))
      string.df$eta.q  [string.df$no.cmpt == 1] <- 0
    if (all(c("no.cmpt", "eta.vp2") %in% names(string.df)))
      string.df$eta.vp2[string.df$no.cmpt < 3] <- 0
    if (all(c("no.cmpt", "eta.q2") %in% names(string.df)))
      string.df$eta.q2 [string.df$no.cmpt < 3] <- 0

    # Michaelis–Menten corrections
    if ("mm" %in% names(string.df)) {
      if (all(c("mm", "eta.vmax") %in% names(string.df)))
        string.df$eta.vmax[string.df$mm == 0] <- 0
      if (all(c("mm", "eta.km")   %in% names(string.df)))
        string.df$eta.km  [string.df$mm == 0] <- 0
      if (all(c("mm", "eta.cl")   %in% names(string.df)))
        string.df$eta.cl  [string.df$mm == 1] <- 0
    }

    # Oral absorption corrections
    oral_cols <-
      c(
        "abs.delay",
        "abs.type",
        "eta.ka",
        "eta.tlag",
        "eta.D2",
        "eta.F1",
        "eta.Fr",
        "eta.mtt",
        "eta.n",
        "eta.bio",
        "abs.bio"
      )
    is_oral <- any(oral_cols %in% names(string.df))
    if (is_oral) {
      if (all(c("abs.delay", "eta.tlag") %in% names(string.df)))
        string.df$eta.tlag[string.df$abs.delay != 1] <- 0
      if ("abs.delay" %in% names(string.df)) {
        if ("eta.mtt" %in% names(string.df))
          string.df$eta.mtt[string.df$abs.delay != 2] <- 0
        if ("eta.n"   %in% names(string.df))
          string.df$eta.n  [string.df$abs.delay != 2] <- 0
        if ("eta.bio" %in% names(string.df))
          string.df$eta.bio[string.df$abs.delay != 2] <- 0
      }
      if (all(c("abs.bio", "abs.type") %in% names(string.df)))
        string.df$abs.bio[string.df$abs.bio == 1 &
                            string.df$abs.type == 4] <- 0
      if (all(c("eta.F1", "abs.bio") %in% names(string.df)))
        string.df$eta.F1[string.df$abs.bio == 0] <- 0
      if (all(c("eta.Fr", "abs.type") %in% names(string.df)))
        string.df$eta.Fr[string.df$abs.type != 4] <- 0
      if (all(c("eta.D2", "abs.type") %in% names(string.df)))
        string.df$eta.D2[string.df$abs.type == 1] <- 0
    }

    # Count IIVs

    # neta1
    cols1 <- intersect(c("eta.vmax", "eta.km"), names(string.df))

    if (length(cols1) == 0) {
      neta1 <- rep(0, nrow(string.df))
    } else {
      neta1 <- rowSums(string.df[, cols1, drop = FALSE], na.rm = TRUE)
    }

    # neta2
    cols2 <- intersect(c("eta.cl", "eta.vc", "eta.vp", "eta.vp2", "eta.q", "eta.q2"),
                       names(string.df))
    if (length(cols2) == 0) {
      neta2 <- rep(0, nrow(string.df))
    } else {
      neta2 <- rowSums(string.df[, cols2, drop = FALSE], na.rm = TRUE)
    }

    # neta3
    cols3 <- intersect(c("eta.ka", "eta.tlag", "eta.D2", "eta.F1", "eta.Fr",
                         "eta.mtt", "eta.n", "eta.bio"),
                       names(string.df))
    if (is_oral && length(cols3) > 0) {
      neta3 <- rowSums(string.df[, cols3, drop = FALSE], na.rm = TRUE)
    } else {
      neta3 <- rep(0, nrow(string.df))
    }


    # mcorr fix
    if ("mcorr" %in% names(string.df)) {
      string.df$mcorr[string.df$mcorr == 1 & neta1 < 2 & neta2 < 2] <- 0
    }

    # Force add IIV if none
    no_iiv <- (neta1 == 0 & neta2 == 0 & neta3 == 0)
    if (any(no_iiv) && "mm" %in% names(string.df)) {
      if ("eta.cl" %in% names(string.df))
        string.df$eta.cl  [no_iiv & string.df$mm == 0] <- 1
      if ("eta.vmax" %in% names(string.df))
        string.df$eta.vmax[no_iiv & string.df$mm == 1] <- 1
    }

    # --- Step 4: If GA, re-encode categorical code back to binary code ---
    if (code.source == "GA") {
      # no.cmpt encoding
      if (string.df$no.cmpt == 1) {
        no.cmpt1 <- 0
        no.cmpt2 <- 1
      } else if (string.df$no.cmpt == 2) {
        no.cmpt1 <- 1
        no.cmpt2 <- 0
      } else if (string.df$no.cmpt == 3) {
        no.cmpt1 <- 1
        no.cmpt2 <- 1
      } else {
        stop("Invalid no.cmpt value.")
      }

      # rv encoding
      if (string.df$rv == 1) {
        rv1 <- 0
        rv2 <- 1
      } else if (string.df$rv == 2) {
        rv1 <- 1
        rv2 <- 0
      } else if (string.df$rv == 3) {
        rv1 <- 1
        rv2 <- 1
      } else {
        stop("Invalid rv value.")
      }

      string.df <- data.frame(
        no.cmpt1 = no.cmpt1,
        no.cmpt2 = no.cmpt2,
        eta.km = string.df$eta.km,
        eta.vc = string.df$eta.vc,
        eta.vp = string.df$eta.vp,
        eta.vp2 = string.df$eta.vp2,
        eta.q = string.df$eta.q,
        eta.q2 = string.df$eta.q2,
        mm = string.df$mm,
        mcorr = string.df$mcorr,
        rv1 = rv1,
        rv2 = rv2
      )

    }

    if (code.source %in% c("ACO", "TS")) {
        drop_cols <- intersect(c("eta.vmax", "eta.cl"), names(string.df))
        string.df <- string.df[, !(names(string.df) %in% drop_cols), drop = FALSE]
    }
    # --- Step 5: Return as vector ---
    return(as.vector(t(string.df)))
  }


#' Parse GA Binary Code to Categorical Code
#'
#' @description
#' Converts Genetic Algorithm (GA) binary model codes into categorical codes
#' that the main pharmacokinetic modeling system can interpret.
#'
#' This function is **only for GA** binary input.
#' - It decodes `no.cmpt1` / `no.cmpt2` into `no.cmpt` (number of compartments).
#' - It decodes `rv1` / `rv2` into `rv` (residual error model type).
#' - Keeps all other parameters as-is.
#' - Supports both `ivcase` (intravenous) and `oralcase` (oral administration).
#'
#' @param string Numeric vector representing GA binary code.
#' @param search.space Character: `"ivcase"` or `"oralcase"`.
#'   Determines the parameter set to be extracted.
#'
#' @return A numeric vector representing the corresponding categorical code.
#' @export
parseCode <- function(string, search.space = "ivbase") {
  # --- Step 1: Ensure numeric vector ---
  string <- as.numeric(as.vector(as.matrix(string)))

  # --- Step 2: Check search space validity ---
  if (!search.space %in% c("ivbase", "oralbase")) {
    stop('search.space must be "ivbase" or "oralbase".')
  }

  # --- Step 3: Extract GA binary parameters ---
  if (search.space == "ivbase") {
    no.cmpt1 <- string[1]
    no.cmpt2 <- string[2]
    eta.km   <- string[3]
    eta.vc   <- string[4]
    eta.vp   <- string[5]
    eta.vp2  <- string[6]
    eta.q    <- string[7]
    eta.q2   <- string[8]
    mm       <- string[9]
    mcorr    <- string[10]
    rv1      <- string[11]
    rv2     <- string[12]
    eta.ka   <- NULL
  } else if (search.space == "oralbase") {
    no.cmpt1 <- string[1]
    no.cmpt2 <- string[2]
    eta.km   <- string[3]
    eta.vc   <- string[4]
    eta.vp   <- string[5]
    eta.vp2  <- string[6]
    eta.q    <- string[7]
    eta.q2   <- string[8]
    eta.ka   <- string[9]
    mm       <- string[10]
    mcorr    <- string[11]
    rv1      <- string[12]
    rv2      <- string[13]
  }

  # --- Step 4: Decode compartment number ---
  if (no.cmpt1 == 0 && no.cmpt2 %in% c(0, 1)) {
    no.cmpt <- 1
  } else if (no.cmpt1 == 1 && no.cmpt2 == 0) {
    no.cmpt <- 2
  } else if (no.cmpt1 == 1 && no.cmpt2 == 1) {
    no.cmpt <- 3
  } else {
    stop("Invalid no.cmpt1/no.cmpt2 combination for GA.")
  }

  # --- Step 5: Decode residual error model type ---
  if (rv1 == 0 && rv2 %in% c(0, 1)) {
    rv <- 1
  } else if (rv1 == 1 && rv2 == 0) {
    rv <- 2
  } else if (rv1 == 1 && rv2 == 1) {
    rv <- 3
  } else {
    stop("Unexpected rv1/rv2 combination for GA.")
  }

  # --- Step 6: Assemble categorical code ---
  string.df <- data.frame(
    no.cmpt = no.cmpt,
    eta.km  = eta.km,
    eta.vc  = eta.vc,
    eta.vp  = eta.vp,
    eta.vp2 = eta.vp2,
    eta.q   = eta.q,
    eta.q2  = eta.q2,
    mm      = mm,
    mcorr   = mcorr,
    rv      = rv,
    stringsAsFactors = FALSE
  )

  # If oral, insert eta.ka after eta.q2
  if (!is.null(eta.ka)) {
    pos <- which(names(string.df) == "eta.q2")
    string.df <- cbind(string.df[, 1:pos, drop = FALSE],
                       eta.ka = eta.ka,
                       string.df[, (pos + 1):ncol(string.df), drop = FALSE])
  }

  # --- Step 7: Return as numeric vector ---
  return(as.vector(t(string.df)))
}
