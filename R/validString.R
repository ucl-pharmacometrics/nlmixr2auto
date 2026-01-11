#' Apply parameter dependency rules
#'
#' Applies dependency constraints among structural and statistical flags in a
#' model-code parameter list to produce a feasible combination.
#'
#' @param params Named list of model-code parameters. Elements are typically
#'   scalar categorical values or 0/1 flags. Unknown elements are ignored.
#'
#' @details
#' Corrections are applied in the following groups:
#' \itemize{
#'   \item Compartment rules: disable peripheral IIV terms when "no.cmpt" implies
#'     they are not used.
#'   \item Michaelis-Menten rules: enable or disable "eta.vmax", "eta.km", and
#'     "eta.cl" based on "mm".
#'   \item Oral absorption rules: enable or disable oral-related terms based on
#'     "abs.delay", "abs.type", and "abs.bio".
#'   \item Correlation rules: disable "mcorr" when too few IIV terms are present.
#'   \item IIV requirement: ensure at least one IIV term is present by enabling a
#'     default term consistent with "mm".
#' }
#'
#' @return A named list with corrected parameter values.
#'
#' @author Zhonghui Huang
#'
#' @examples
#' params <- list(
#'   no.cmpt = 1, mm = 0, mcorr = 1,
#'   eta.vc = 1, eta.cl = 0, eta.vp = 1, eta.q = 1
#' )
#' applyParamDeps(params)
#'
#' params2 <- list(
#'   no.cmpt = 2, mm = 1,
#'   eta.vmax = 0, eta.km = 0, eta.cl = 1
#' )
#' applyParamDeps(params2)
#'
#' @export

applyParamDeps <- function(params) {
  string.df <-
    as.data.frame(as.list(params), stringsAsFactors = FALSE)
  # Compartment-based corrections
  if (all(c("no.cmpt", "eta.vp") %in% names(string.df))) {
    string.df$eta.vp[string.df$no.cmpt == 1] <- 0
  }
  if (all(c("no.cmpt", "eta.q") %in% names(string.df))) {
    string.df$eta.q[string.df$no.cmpt == 1] <- 0
  }
  if (all(c("no.cmpt", "eta.vp2") %in% names(string.df))) {
    string.df$eta.vp2[string.df$no.cmpt < 3] <- 0
  }
  if (all(c("no.cmpt", "eta.q2") %in% names(string.df))) {
    string.df$eta.q2[string.df$no.cmpt < 3] <- 0
  }

  # Michaelis-Menten corrections
  if ("mm" %in% names(string.df)) {
    if (all(c("mm", "eta.vmax") %in% names(string.df))) {
      string.df$eta.vmax[string.df$mm == 0] <- 0
    }
    if (all(c("mm", "eta.km") %in% names(string.df))) {
      string.df$eta.km[string.df$mm == 0] <- 0
    }
    if (all(c("mm", "eta.cl") %in% names(string.df))) {
      string.df$eta.cl[string.df$mm == 1] <- 0
    }
  }

  # Oral absorption corrections
  oral_cols <- c(
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
  has_oral_params <- any(oral_cols %in% names(string.df))

  if (has_oral_params) {
    if (all(c("abs.delay", "eta.tlag") %in% names(string.df))) {
      string.df$eta.tlag[string.df$abs.delay != 1] <- 0
    }
    if ("abs.delay" %in% names(string.df)) {
      if ("eta.mtt" %in% names(string.df)) {
        string.df$eta.mtt[string.df$abs.delay != 2] <- 0
      }
      if ("eta.n" %in% names(string.df)) {
        string.df$eta.n[string.df$abs.delay != 2] <- 0
      }
      if ("eta.bio" %in% names(string.df)) {
        string.df$eta.bio[string.df$abs.delay != 2] <- 0
      }
    }
    if (all(c("abs.bio", "abs.type") %in% names(string.df))) {
      string.df$abs.bio[string.df$abs.bio == 1 &
                          string.df$abs.type == 4] <- 0
    }
    if (all(c("eta.F1", "abs.bio") %in% names(string.df))) {
      string.df$eta.F1[string.df$abs.bio == 0] <- 0
    }
    if (all(c("eta.Fr", "abs.type") %in% names(string.df))) {
      string.df$eta.Fr[string.df$abs.type != 4] <- 0
    }
    if (all(c("eta.D2", "abs.type") %in% names(string.df))) {
      string.df$eta.D2[string.df$abs.type == 1] <- 0
    }
  }

  # Count IIVs
  cols1 <- intersect(c("eta.vmax", "eta.km"), names(string.df))
  neta1 <- if (length(cols1) == 0) {
    rep(0, nrow(string.df))
  } else {
    rowSums(string.df[, cols1, drop = FALSE], na.rm = TRUE)
  }

  cols2 <-
    intersect(c("eta.cl", "eta.vc", "eta.vp", "eta.vp2", "eta.q", "eta.q2"),
              names(string.df))
  neta2 <- if (length(cols2) == 0) {
    rep(0, nrow(string.df))
  } else {
    rowSums(string.df[, cols2, drop = FALSE], na.rm = TRUE)
  }

  cols3 <- intersect(
    c(
      "eta.ka",
      "eta.tlag",
      "eta.D2",
      "eta.F1",
      "eta.Fr",
      "eta.mtt",
      "eta.n",
      "eta.bio"
    ),
    names(string.df)
  )
  neta3 <- if (has_oral_params && length(cols3) > 0) {
    rowSums(string.df[, cols3, drop = FALSE], na.rm = TRUE)
  } else {
    rep(0, nrow(string.df))
  }

  # Correlation correction
  if ("mcorr" %in% names(string.df)) {
    string.df$mcorr[string.df$mcorr == 1 & neta1 < 2 & neta2 < 2] <- 0
  }

  # Force at least one IIV
  no_iiv <- (neta1 == 0 & neta2 == 0 & neta3 == 0)
  if (any(no_iiv) && "mm" %in% names(string.df)) {
    if ("eta.cl" %in% names(string.df)) {
      string.df$eta.cl[no_iiv & string.df$mm == 0] <- 1
    }
    if ("eta.vmax" %in% names(string.df)) {
      string.df$eta.vmax[no_iiv & string.df$mm == 1] <- 1
    }
  }
  return(as.list(string.df[1, ]))
}



#' Validate and correct model string for ACO/TS
#'
#' Validates model parameter strings from ACO or tabu search algorithms.
#'
#' @param string Numeric vector representing categorical model encoding.
#' @param search.space Character string specifying which search space to use.
#'   Options are "ivbase", "oralbase", or "custom". Default is "ivbase".
#' @param custom_config List, configuration for custom search spaces. Required
#'   when search.space is "custom".
#'
#' @details
#' The input string is interpreted using the parameter order defined by the
#' selected search space configuration (for "custom", this is custom_config$params).
#' The function:
#' \enumerate{
#'   \item Maps the input vector to named parameters via parseParams().
#'   \item Enforces model constraints via applyParamDeps().
#'   \item Returns only the parameters that belong to the search space, using
#'     the same order as space_cfg$params.
#' }
#' This design ensures the returned vector is compatible with downstream model
#' generation and with binary encoding wrappers (for example, validStringbinary()).
#'
#' @return Numeric vector of validated and corrected parameters (categorical).
#'
#' @author Zhonghui Huang
#'
#' @examples
#' # Example 1: ivbase, 1 compartment disables peripheral terms.
#' invalid_iv <- c(1, 1, 1, 1, 0, 1, 0, 0, 0, 1)
#' validStringcat(invalid_iv, "ivbase")
#'
#' # Example 2: oralbase, mm = 0 forces eta.km to 0.
#' invalid_oral <- c(2, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1)
#' validStringcat(invalid_oral, "oralbase")
#'
#' # Example 3: custom, mcorr is cleared when there are not enough IIV terms.
#' simple_config <- list(
#'   route = "bolus",
#'   params = c("eta.vc", "mcorr", "rv"),
#'   param_dependencies = list(),
#'   fixed_params = list(no.cmpt = 1, eta.cl = 1, allometric_scaling = 1)
#' )
#' invalid_custom <- c(0, 1, 4)
#' validStringcat(invalid_custom, "custom", custom_config = simple_config)
#'
#' @seealso
#' \code{\link{validStringbinary}} for the GA wrapper using binary encoding.
#' \code{\link{parseParams}} for mapping vectors to named parameters.
#' \code{\link{applyParamDeps}} for constraint enforcement rules.
#'
#' @export
validStringcat <- function(string,
                           search.space = "ivbase",
                           custom_config = NULL) {
  if (missing(custom_config)) {
    custom_config <- NULL
  }

  if (search.space == "ivbase") {
    space_cfg <- spaceConfig("ivbase")
  } else if (search.space == "oralbase") {
    space_cfg <- spaceConfig("oralbase")
  } else if (search.space == "custom") {
    if (is.null(custom_config)) {
      stop("custom_config must be provided when search.space is 'custom'.",
           call. = FALSE)
    }
    space_cfg <- custom_config
  } else {
    stop("search.space must be 'ivbase', 'oralbase', or 'custom'.",
         call. = FALSE)
  }
  if (is.null(space_cfg$params) || length(space_cfg$params) == 0) {
    stop("Configuration must provide a non-empty 'params' vector.",
         call. = FALSE)
  }
  # Parse parameters and apply dependency rules
  params <- parseParams(string, space_cfg)
  params <- applyParamDeps(params)

  # Return corrected categorical string
  if (search.space == "ivbase") {
    out <- unlist(params[space_cfg$params], use.names = FALSE)
  } else if (search.space == "oralbase") {
    out <- unlist(params[space_cfg$params], use.names = FALSE)
  } else {
    # search.space == "custom"
    out <- unlist(params[space_cfg$params], use.names = FALSE)
  }
  if (any(is.na(out))) {
    missing_names <- space_cfg$params[is.na(out)]
    stop(sprintf(
      "Missing parameters after correction: %s",
      paste(missing_names, collapse = ", ")
    ),
    call. = FALSE)
  }
  return(as.numeric(out))
}


#' Validate and correct model string for GA
#'
#' Validates model parameter strings from genetic algorithms.
#'
#' @param string Numeric vector representing binary model encoding (0/1).
#' @param search.space Character string specifying which search space to use.
#'   Options are "ivbase", "oralbase", or "custom". Default is "ivbase".
#' @param custom_config List, configuration for custom search spaces. Required
#'   when search.space is "custom".
#'
#' @details
#' The input string is a binary chromosome (0/1). The function:
#' \enumerate{
#'   \item Decodes the binary chromosome to a categorical parameter vector using
#'     \code{decodeBinary}.
#'   \item Applies model constraints in categorical space by calling
#'     \code{validStringcat}.
#'   \item Encodes the corrected categorical vector back to binary using
#'     \code{encodeBinary}.
#' }
#' This design keeps all correction rules in \code{validStringcat} and
#' makes the GA version a thin wrapper around the categorical validator.
#'
#' @return Numeric vector of validated and corrected parameters (binary encoding).
#'
#' @author Zhonghui Huang
#'
#' @examples
#' # Example 1: ivbase, 1 compartment disables peripheral terms.
#' # Bits 1-2 encode no.cmpt; here 00 maps to 1.
#' invalid_iv <- c(0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1)
#' validStringbinary(invalid_iv, "ivbase")
#'
#' # Example 2: oralbase, mm = 0 forces eta.km to 0.
#' # Bits 12-13 encode rv for oralbase.
#' invalid_oral <- c(1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1)
#' validStringbinary(invalid_oral, "oralbase")
#'
#' # Example 3: custom, mcorr is cleared when there are not enough IIV terms.
#' simple_config <- list(
#'   route = "bolus",
#'   params = c("eta.vc", "mcorr", "rv"),
#'   param_dependencies = list(),
#'   fixed_params = list(no.cmpt = 1, eta.cl = 1, allometric_scaling = 1)
#' )
#' # custom encoding: eta.vc (1 bit), mcorr (1 bit), rv (2 bits)
#' invalid_custom <- c(0, 1, 1, 1)  # eta.vc=0, mcorr=1, rv=4
#' validStringbinary(invalid_custom, "custom", custom_config = simple_config)
#'
#' @seealso
#' \code{\link{validStringcat}} for categorical validation used by ACO/TS.
#' \code{\link{decodeBinary}} and \code{\link{encodeBinary}} for encoding
#' conversions.
#'
#' @export

validStringbinary <- function(string,
                              search.space = "ivbase",
                              custom_config = NULL) {
  if (missing(custom_config)) {
    custom_config <- NULL
  }
  # Decode binary -> categorical
  categorical <- decodeBinary(
    binary_string = string,
    search.space = search.space,
    custom_config = custom_config
  )

  # Validate/correct in categorical space (reuse ACO/TS logic)
  if (search.space == "ivbase") {
    corrected_categorical <- validStringcat(
      string = categorical,
      search.space = "ivbase",
      custom_config = custom_config
    )
  } else if (search.space == "oralbase") {
    corrected_categorical <- validStringcat(
      string = categorical,
      search.space = "oralbase",
      custom_config = custom_config
    )
  } else if (search.space == "custom") {
    corrected_categorical <- validStringcat(
      string = categorical,
      search.space = "custom",
      custom_config = custom_config
    )
  } else {
    stop("search.space must be 'ivbase', 'oralbase', or 'custom'.",
         call. = FALSE)
  }

  # Encode categorical -> binary
  binary_out <- encodeBinary(
    categorical_string = corrected_categorical,
    search.space = search.space,
    custom_config = custom_config
  )

  return(as.numeric(binary_out))
}
