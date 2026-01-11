#' Parse model coding vector to model name
#'
#' Converts an ordinal model-coding vector into a single pharmacokinetic model
#' name string, using a search space configuration.
#'
#' @param modcode Numeric vector of model-coding flags/values in the order
#'   defined by the search space configuration.
#' @param search.space Character string specifying which search space to use.
#'   Options are "ivbase", "oralbase", or "custom". Default is "ivbase".
#' @param custom_config Optional named list defining a custom parameter structure.
#'   If provided, the parameter names are taken from the names of this list.
#'   If NULL, a default parameter structure is used based on the selected
#'   search space.
#' @details
#' The function selects a configuration (either a custom configuration when
#' search.space is "custom", or the predefined configuration from
#' spaceConfig(search.space) otherwise). It then decodes modcode with
#' parseParams() and assembles an underscore-separated model name.
#'
#' The name is built from these blocks in order: prefix, compartments,
#' optional absorption, ETA block, elimination, correlation, residual error,
#' and optional allometric scaling.
#'
#' Key rules:
#' \itemize{
#'   \item The prefix is taken from config$prefix when available; otherwise from
#'         config$route; otherwise from search.space.
#'   \item The absorption block is included when any absorption option is present.
#'         If abs.type, abs.delay, and abs.bio are all missing/NA, the absorption
#'         block is included only for oral or mixed routes using FO_abs; otherwise
#'         it is omitted.
#'   \item The ETA block includes all eta.* terms equal to 1 (NA values are ignored),
#'         and also forces Vmax when mm equals 1; otherwise it forces CL.
#'   \item Elimination uses FO_elim when mm is 0 or NA, and MM_elim when mm is 1.
#'         Correlation uses uncorrelated when mcorr is 0 or NA, and correlated when
#'         mcorr is 1. Residual error uses add, prop, or combined.
#'   \item Allometric scaling is omitted when allometric_scaling is 0 or NA; otherwise
#'         it appends "asWT", "asBMI", or "asFFM".
#' }
#'
#' @return A character string representing the constructed model name.
#'
#' @author Zhonghui Huang
#'
#' @examples
#' # Example 1: Parse IV base model name
#' parseName(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 1), "ivbase")
#'
#' # Example 2: Parse oral base model name
#' parseName(c(2, 1, 1, 0, 0, 1, 1, 1, 0, 1, 3), "oralbase")
#'
#' # Example 3: Parse custom configuration model name
#' custom_config <- list(
#'   prefix = "custom",
#'   route  = "oral",
#'   params = c("no.cmpt", "eta.cl", "eta.vc", "mm", "mcorr", "rv"),
#'   param_dependencies = list(),
#'   fixed_params = list()
#' )
#' parseName(c(2, 1, 0, 0, 1, 2), search.space = "custom",
#' custom_config = custom_config)
#'
#' @seealso \code{spaceConfig()}, \code{parseParams()}
#' @export
#'
parseName <- function(modcode,
                      search.space = NULL,
                      custom_config = NULL) {

  if (missing(custom_config)){
    custom_config <-NULL
  }
  if (search.space == "custom" && !is.null(custom_config)) {
    config <- custom_config
  } else {
    config <- spaceConfig(search.space)
  }

  # ---- parse modcode -> named params ----
  params <- parseParams(string = modcode,config = config)

  # ---- prefix ----
  prefix <- if (!is.null(config$prefix)) {
    config$prefix
  } else if (!is.null(config$route)) {
    config$route
  } else if (!is.null(search.space)) {
    as.character(search.space)
  } else {
    "model"
  }

  # ---- extract features ----
  no_cmpt <-
    if ("no.cmpt" %in% names(params))
      as.numeric(params[["no.cmpt"]])
  else
    NA_real_
  cmpt_string <-
    if (!is.na(no_cmpt))
      paste0(no_cmpt, "cmpt")
  else
    "cmptNA"

  mm    <-
    if ("mm" %in% names(params))
      as.numeric(params[["mm"]])
  else
    NA_real_
  mcorr <-
    if ("mcorr" %in% names(params))
      as.numeric(params[["mcorr"]])
  else
    NA_real_
  rv    <-
    if ("rv" %in% names(params))
      as.numeric(params[["rv"]])
  else
    NA_real_

  # ---- absorption features ----
  abs_type  <-
    if ("abs.type" %in% names(params))
      as.numeric(params[["abs.type"]])
  else
    NA_real_
  abs_delay <-
    if ("abs.delay" %in% names(params))
      as.numeric(params[["abs.delay"]])
  else
    NA_real_
  abs_bio   <-
    if ("abs.bio" %in% names(params))
      as.numeric(params[["abs.bio"]])
  else
    NA_real_

  route_val <-
    if ("route" %in% names(params)) {
      as.character(params[["route"]])
    } else if (!is.null(config$route)) {
      as.character(config$route)
    } else {
      NA_character_
    }

  abs_type_string <- if (is.na(abs_type)) {
    "FO"
  } else {
    switch(
      as.character(abs_type),
      "1" = "FOabs",
      "2" = "ZOabs",
      "3" = "SEQabs",
      "4" = "MIXabs",
      paste0("T", abs_type)
    )
  }

  abs_delay_string <- if (is.na(abs_delay)) {
    NA_character_
  } else {
    switch(
      as.character(abs_delay),
      "0" = "d0",
      "1" = "TLAG",
      "2" = "TR",
      paste0("D", abs_delay)
    )
  }

  abs_bio_string <- if (is.na(abs_bio)) {
    NA_character_
  } else {
    switch(as.character(abs_bio),
           "0" = "F0",
           "1" = "F1",
           paste0("F", abs_bio))
  }

  abs_string <-
    if (is.na(abs_type) && is.na(abs_delay) && is.na(abs_bio)) {
      if (!is.na(route_val) &&
          route_val %in% c("oral", "mixed_iv_oral")) {
        paste(c("FO", "abs"), collapse = "_")
      } else {
        NULL
      }
    } else {
      abs_pieces <- c(abs_type_string, "abs")
      if (!is.na(abs_delay_string))
        abs_pieces <- c(abs_pieces, abs_delay_string)
      if (!is.na(abs_bio_string))
        abs_pieces <- c(abs_pieces, abs_bio_string)
      paste(abs_pieces, collapse = "_")
    }

  eta_names <- grep("^eta\\.", names(params), value = TRUE)
  eta_vals <-
    suppressWarnings(as.numeric(unlist(params[eta_names])))
  eta_active <- eta_names[!is.na(eta_vals) & eta_vals == 1]

  key <- sub("^eta\\.", "", eta_active)
  key_l <- tolower(key)

  map <- c(
    "cl" = "CL",
    "vc" = "VC",
    "vp" = "VP",
    "vp2" = "VP2",
    "q"  = "Q",
    "q2" = "Q2",
    "km" = "KM",
    "vmax" = "Vmax",
    "ka" = "KA",
    "tlag" = "TLAG",
    "n" = "N",
    "mtt" = "MTT",
    "bio" = "BIO",
    "d2" = "D2",
    "f1" = "F1",
    "fr" = "FR"
  )

  mapped <- unname(map[key_l])
  eta_labels <- ifelse(is.na(mapped), toupper(key_l), mapped)

  forced <- if (!is.na(mm) && mm == 1)
    "Vmax"
  else
    "CL"
  eta_parts <- unique(c(forced, eta_labels))
  eta_string <- paste0("eta", paste(eta_parts, collapse = ""))

  mm_string <- if (is.na(mm)) {
    "FOelim"
  } else if (mm == 0) {
    "FOelim"
  } else {
    "MMelim"
  }

  corr_string <-
    if (is.na(mcorr) || mcorr == 0)
      "uncorrelated"
  else
    "correlated"

  rv_string <- if (is.na(rv)) {
    "add"
  } else {
    switch(
      as.character(rv),
      "1" = "add",
      "2" = "prop",
      "3" = "combined",
      paste0("rv", rv)
    )
  }

  allo <-
    if ("allometric_scaling" %in% names(params))
      as.numeric(params[["allometric_scaling"]])
  else
    NA_real_

  allo_string <- if (is.na(allo) || allo == 0) {
    NULL
  } else {
    switch(
      as.character(allo),
      "1" = "asWT",
      "2" = "asBMI",
      "3" = "asFFM",
      paste0("allo", allo)
    )
  }

  parts <- list(
    prefix,
    cmpt_string,
    abs_string,
    eta_string,
    mm_string,
    corr_string,
    rv_string,
    allo_string
  )

  parts <- parts[!vapply(parts, is.null, logical(1))]
  paste(unlist(parts), collapse = "_")
}
