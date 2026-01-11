#' Parse string vector to model parameters
#'
#' Converts a numeric vector of parameter values into a named list of model
#' parameters based on the search space configuration.
#'
#' @param string Numeric vector containing parameter values in the order
#'   specified by the search space configuration
#' @param config List object returned by \code{spaceConfig()},
#'   containing parameter definitions and dependencies
#'
#' @details
#' This function performs three main operations:
#' \enumerate{
#'   \item Maps the input vector to named parameters
#'   \item Computes dependent parameter values using defined functions
#'   \item Adds fixed parameters and route information
#' }
#'
#' @return A named list containing:
#'   \itemize{
#'     \item All parameters specified in config$params with their values
#'     \item Computed dependent parameters based on param_dependencies
#'     \item Fixed parameters from `fixed_params`
#'     \item Administration route from config$route
#'   }
#'
#' @author Zhonghui Huang
#'
#' @examples
#' # Example 1: Parse IV base model parameters
#' config_iv <- spaceConfig("ivbase")
#' parseParams(c(2, 1, 1, 0, 0, 1, 1, 0, 1, 1), config_iv)
#'
#' # Example 2: Parse oral base model parameters
#' config_oral <- spaceConfig("oralbase")
#' parseParams(c(2, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1), config_oral)
#'
#' # Example 3: Parse custom configuration parameters
#' custom_config <- list(route = "oral", params = c("no.cmpt", "eta.cl", "eta.vc"),
#'                       param_dependencies = list(), fixed_params = list(mm = 0))
#' parseParams(c(1, 1, 1), custom_config)
#'
#' @seealso \code{spaceConfig()}, \code{mod.run()}
#' @export
#'
parseParams <- function(string, config) {
  # Validate custom search space
  if (is.null(config$params)) {
    stop("Custom search space requires explicit parameter specification",
         call. = FALSE)
  }

  if (length(string) != length(config$params)) {
    stop(
      sprintf(
        "Length mismatch: string has %d elements but search space expects %d parameters (%s)",
        length(string),
        length(config$params),
        paste(config$params, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  params <- stats::setNames(as.list(string), config$params)

  if (length(config$param_dependencies) > 0) {
    for (dep_name in names(config$param_dependencies)) {
      dep_func <- config$param_dependencies[[dep_name]]
      # Extract arguments needed by the dependency function
      dep_args <- lapply(names(formals(dep_func)), function(x) {
        as.numeric(params[[x]])
      })
      params[[dep_name]] <- do.call(dep_func, dep_args)
    }
  }

  if (length(config$fixed_params) > 0) {
    params <- c(params, config$fixed_params)
  }

  if (!is.null(config$route)) {
    params$route <- config$route
  }

  # Define standard parameter order
  standard_order <- c(
    "no.cmpt",
    "eta.vmax",
    "eta.km",
    "eta.cl",
    "eta.vc",
    "eta.vp",
    "eta.vp2",
    "eta.q",
    "eta.q2",
    "eta.ka",
    "eta.tlag",
    "eta.n",
    "eta.mtt",
    "eta.bio",
    "eta.D2",
    "eta.F1",
    "eta.Fr",
    "mm",
    "mcorr",
    "rv",
    "allometric_scaling",
    "abs.type",
    "abs.delay",
    "abs.bio",
    "route"
  )

  param_names <- names(params)
  in_standard <- param_names[param_names %in% standard_order]
  not_in_standard <- param_names[!param_names %in% standard_order]
  ordered_standard <- standard_order[standard_order %in% in_standard]
  final_order <- c(ordered_standard, sort(not_in_standard))
  params <- params[final_order]

  return(params)
}
