#' Create a base model code for single-start model search algorithms
#'
#' Constructs a named numeric vector defining the initial structural and
#' inter-individual variability model configuration used in single-start
#' automated PK model search algorithms.
#'
#' @param search.space Character, one of "ivbase" or "oralbase". Default is "ivbase".
#' @details
#' Two search spaces are supported: "ivbase" and "oralbase". A user-specified
#' initial model code can be provided via the custom_base argument. The input is
#' validated for numerical type and expected length, and standardized element
#' names are applied before returning. The function is currently used in
#' stepwise selection and tabu search routines, where a single starting model
#' is iteratively updated.
#'
#' @return
#' For search.space = "ivbase": a named integer vector of length 9 containing:
#' - no.cmpt   - Number of compartments
#' - eta.km    - IIV flag for \eqn{K_m}
#' - eta.vc    - IIV flag for \eqn{V_c}
#' - eta.vp    - IIV flag for \eqn{V_p}
#' - eta.vp2   - IIV flag for \eqn{V_{p2}}
#' - eta.q     - IIV flag for \eqn{Q}
#' - eta.q2    - IIV flag for \eqn{Q_2}
#' - mm        - Michaelis–Menten term flag
#' - mcorr     - Correlation flag among ETAs
#' - rv        - Residual error model code
#'
#' For search.space = "oralbase": a named integer vector of length 11, including
#' all fields above plus:
#' - eta.ka    - IIV flag for \eqn{k_a} (oral absorption rate constant)
#'
#' @author Zhonghui Huang
#'
#' @examples
#' base_model("ivbase")
#' base_model("oralbase")
#' @export
base_model <- function(search.space = "ivbase") {
  # Element names for the full oralbase vector
  element_names <- c(
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
  )

  default_code <- c(
    no.cmpt = 1,
    eta.km  = 0,
    eta.vc  = 0,
    eta.vp  = 0,
    eta.vp2 = 0,
    eta.q   = 0,
    eta.q2  = 0,
    eta.ka  = 0,
    mm      = 0,
    mcorr   = 0,
    rv      = 3
  )

  if (search.space == "ivbase") {
    return(default_code[c(1:7, 9:11)])
  } else if (search.space == "oralbase") {
    return(default_code)
  } else {
    stop("Invalid search.space. Must be 'ivbase' or 'oralbase'.")
  }
}
