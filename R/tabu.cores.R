#' Generate neighbor models
#'
#' Generates a set of neighbor models from a given current model code.
#' The neighborhood is defined as all single-variable changes
#' (1-bit or 1-step modifications). Each generated neighbor is passed
#' through \code{validateModels()} to ensure the resulting code
#' corresponds to a valid pharmacometric model.
#'
#' Optionally, the function can restrict the number of neighbors
#' by random sampling (candidate list strategy).
#'
#' @param current_string A named numeric vector representing the current
#'   model code. Names correspond to model features (e.g. \code{no.cmpt},
#'   \code{eta.vc}, \code{rv}), and values to their current states.
#' @param search.space Character scalar: either \code{"ivbase"} or
#'   \code{"oralbase"}. Determines which model features are available.
#' @param candidate.size Integer (optional). Maximum number of neighbors
#'   to return. If \code{NULL} (default), the full neighborhood is returned.
#'   If specified, a random subset of this size is sampled.
#' @param seed Optional random seed for reproducibility when sampling neighbors.
#'
#' @return A \code{data.frame} of valid neighbor model codes, with
#'   the same column names as \code{current_string}.
#'
#' @examples
#' current_string <- c(no.cmpt = 2, eta.km = 0, eta.vc = 1,
#'                     eta.vp = 0, eta.vp2 = 0, eta.q = 1,
#'                     eta.q2 = 0, mm = 0, mcorr = 1, rv = 2)
#' neighbors <- generate_neighbors_df(current_string, search.space = "ivbase")
#' head(neighbors)
#'
#' @export
generate_neighbors_df <- function(current_string,
                                  search.space = c("ivbase", "oralbase"),
                                  candidate.size = NULL,
                                  seed = NULL) {
  search.space <- match.arg(search.space)

  # Define available option sets depending on search space
  if (search.space == "ivbase") {
    options_list <- list(
      no.cmpt = 1:3,
      eta.km  = 0:1,
      eta.vc  = 0:1,
      eta.vp  = 0:1,
      eta.vp2 = 0:1,
      eta.q   = 0:1,
      eta.q2  = 0:1,
      mm      = 0:1,
      mcorr   = 0:1,
      rv      = 1:3
    )
  } else {
    options_list <- list(
      no.cmpt = 1:3,
      eta.km  = 0:1,
      eta.vc  = 0:1,
      eta.vp  = 0:1,
      eta.vp2 = 0:1,
      eta.q   = 0:1,
      eta.q2  = 0:1,
      eta.ka  = 0:1,
      mm      = 0:1,
      mcorr   = 0:1,
      rv      = 1:3
    )
  }

  # Initialize empty neighbor table
  neighbors <-
    data.frame(matrix(ncol = length(current_string), nrow = 0))
  colnames(neighbors) <- names(current_string)

  # Generate neighbors by modifying one variable at a time
  for (variable in names(current_string)) {
    options <- options_list[[variable]]
    for (option in options) {
      if (option != current_string[variable]) {
        new_neighbor <- current_string
        new_neighbor[variable] <- option
        neighbors <- rbind(neighbors, new_neighbor)
      }
    }
  }

  # Validate neighbors: ensures only legal pharmacometric models remain
  neighbors <- t(vapply(seq_len(nrow(neighbors)),
                        function(i) {
                          validateModels(
                            string = as.numeric(neighbors[i,]),
                            search.space = search.space,
                            code.source = "TS"  # flag to distinguish Tabu Search
                          )
                        },
                        numeric(length(current_string))))

  neighbors <- as.data.frame(neighbors)
  colnames(neighbors) <- names(current_string)

  # Deduplicate after validation
  neighbors <- neighbors[!duplicated(neighbors), , drop = FALSE]

  # Candidate list sampling (optional)
  if (!is.null(candidate.size) &&
      nrow(neighbors) > candidate.size) {
    if (!is.null(seed))
      set.seed(seed)
    neighbors <-
      neighbors[sample(1:nrow(neighbors), candidate.size), , drop = FALSE]
  }

  return(neighbors)
}


#' Detect the move between two solutions
#'
#' @description
#' Given a current solution vector and a neighbor solution vector,
#' this function identifies the specific move (bit change) required to
#' transform the current solution into the neighbor solution.
#'
#' @param prev_string Numeric vector representing the current solution
#'   (e.g., the local best).
#' @param new_string Numeric vector representing a candidate neighbor solution.
#'
#' @return A named list with elements:
#'   \itemize{
#'     \item \code{element} — The variable name that was changed.
#'     \item \code{from} — The original value in the current solution.
#'     \item \code{to} — The new value in the neighbor solution.
#'   }
#'
#' @examples
#' prev <- c(no.cmpt = 2, eta.vc = 0, rv = 2)
#' new  <- c(no.cmpt = 3, eta.vc = 0, rv = 2)
#' detect_move(prev, new)
#'
#' @export
detect_move <- function(prev_string, new_string) {
  diff_idx <- which(prev_string != new_string)
  if (length(diff_idx) == 0) {
    return(NULL)  # no change
  }
  if (length(diff_idx) > 1) {
    warning("Multiple changes detected; only the first will be reported.")
    diff_idx <- diff_idx[1]
  }
  element <- names(prev_string)[diff_idx]
  from <- prev_string[diff_idx]
  to <- new_string[diff_idx]
  return(list(element = element, from = from, to = to))
}


#' Check if a move is tabu
#'
#' @description
#' Given a move (variable, from-value, to-value) and a tabu list,
#' this function checks whether the move is currently forbidden
#' by the tabu list.
#'
#' @param move A list as returned by \code{\link{detect_move}},
#'   containing \code{element}, \code{from}, and \code{to}.
#' @param tabu_list Data frame of tabu elements, with columns:
#'   \code{elements} (variable name), \code{elements.value} (forbidden value),
#'   and \code{tabu.iteration.left} (remaining tabu tenure).
#'
#' @return Logical scalar: \code{TRUE} if the move is tabu, \code{FALSE} otherwise.
#'
#' @examples
#' move <- list(element = "no.cmpt", from = 2, to = 3)
#' tabu_list <- data.frame(
#'   elements = c("no.cmpt", "eta.vc"),
#'   elements.value = c(3, 1),
#'   tabu.iteration.left = c(2, 1)
#' )
#' is_move_tabu(move, tabu_list)
#'
#' @export
is_move_tabu <- function(move, tabu_list) {
  if (is.null(move) || is.null(tabu_list) || nrow(tabu_list) == 0) {
    return(FALSE)
  }
  return(any(
    tabu_list$elements == move$element &
      tabu_list$elements.value == move$to &
      tabu_list$tabu.iteration.left > 0
  ))
}

#' Apply 1-bit perturbation to escape local optimum
#'
#' This function randomly flips one parameter in the current model string
#' to generate a new valid model. Used when the search stagnates.
#'
#' @param prev_string A named numeric vector representing the current model.
#' @param search.space The search space ("ivbase" or "oralbase").
#' @param max.try Maximum number of attempts to find a valid perturbed model.
#'
#' @return A perturbed model string (named numeric vector).
#' @export
perturb_1bit <- function(prev_string, search.space, max.try = 1000) {
  for (i in 1:max.try) {
    disturbed <- prev_string
    flip_idx <- sample(seq_along(prev_string), 1)

    param <- names(prev_string)[flip_idx]
    if (param %in% c("no.cmpt", "rv")) {
      all_options <- 1:3
    } else {
      all_options <- 0:1
    }

    current_val <- as.numeric(disturbed[flip_idx])
    candidates <- setdiff(all_options, current_val)

    if (length(candidates) > 0) {
      disturbed[flip_idx] <- sample(candidates, 1)
    }

    disturbed <- validateModels(string       = disturbed,
                                search.space = search.space,
                                code.source  = "TS")
    names(disturbed) <- names(prev_string)

    if (!all(disturbed == prev_string)) {
      return(disturbed)
    }
  }
  warning("1-bit perturbation failed after ",
          max.try,
          " attempts, returning original string.")
  return(prev_string)
}
