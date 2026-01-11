#' Get search space configuration
#'
#' Retrieve the configuration for a specified search space.
#'
#' @param search.space Character, one of "ivbase" or "oralbase".
#'   Default is "ivbase".
#'
#' @details
#' Pre-defined search spaces:
#' \itemize{
#'   \item "ivbase": IV bolus model, 11 parameters, supports 1 to 3 compartments.
#'   \item "oralbase": Oral model, 12 parameters (adds eta.ka), supports 1 to 3 compartments.
#' }
#'
#' For "ivbase" and "oralbase", param_dependencies handle the relationship between
#' Michaelis-Menten elimination (mm) and the associated variability parameters
#' (eta.vmax, eta.cl).
#' @return A list with four elements:
#' \itemize{
#'   \item route: Administration route ("bolus", "oral", or NULL).
#'   \item params: Character vector of parameter names expected in the string vector.
#'   \item param_dependencies: Named list of functions that compute dependent parameters.
#'   \item fixed_params: Named list of fixed parameter values.
#' }
#'
#' @author Zhonghui Huang
#'
#' @examples
#' # Get IV base configuration
#' config <- spaceConfig("ivbase")
#' config$params
#'
#' # Get oral base configuration
#' config <- spaceConfig("oralbase")
#' config$params
#'
#' @seealso
#' \link{mod.run} for the main function that uses these configurations.
#' \link{parseParams} for parameter parsing using configurations.
#'
#' @export
spaceConfig <- function(search.space = c("ivbase", "oralbase")) {
  search.space <- match.arg(search.space)
  configs <- list(
    ivbase = list(
      route = "bolus",
      params = c(
        "no.cmpt",
        "eta.km",
        "eta.vc",
        "eta.vp",
        "eta.vp2",
        "eta.q",
        "eta.q2",
        "mm",
        "mcorr",
        "rv"
      ),
      param_dependencies = list(
        eta.vmax = function(mm)
          if (mm == 0)
            0
        else
          1,
        eta.cl = function(mm)
          if (mm == 1)
            0
        else
          1
      ),
      # Fixed parameter values
      fixed_params = list()
    ),

    oralbase = list(
      route = "oral",
      params = c(
        "no.cmpt",
        "eta.km",
        "eta.vc",
        "eta.vp",
        "eta.vp2",
        "eta.q",
        "eta.q2",
        "eta.ka",
        "mm",
        "mcorr",
        "rv"
      ),
      param_dependencies = list(
        eta.vmax = function(mm)
          if (mm == 0)
            0
        else
          1,
        eta.cl = function(mm)
          if (mm == 1)
            0
        else
          1
      ),
      fixed_params = list()
    ),

    custom = list(
      route = NULL,
      params = NULL,
      param_dependencies = list(),
      fixed_params = list()
    )
  )

  if (!search.space %in% names(configs)) {
    stop(
      sprintf(
        "Unknown search.space: '%s'. Available options: %s",
        search.space,
        paste(names(configs), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  return(configs[[search.space]])
}
