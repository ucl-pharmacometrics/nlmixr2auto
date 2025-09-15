#' Generate neighbor models
#'
#' Generates a set of neighbor models from a given current model code.
#' The neighborhood is defined as all single-variable changes
#' (1-bit or 1-step modifications). Each generated neighbor is passed
#' through \code{validateModels()} to ensure the resulting code
#' corresponds to a valid pharmacometric model.
#'
#' For each neighbor, both the \emph{original} (pre-validation) and the
#' \emph{validated} (post-validation) codes are retained. This allows
#' downstream functions (e.g. \code{detect_move()}) to distinguish
#' between the intended primary modification and any secondary
#' adjustments introduced by validation.
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
#' @return A \code{list} with two components:
#'   \describe{
#'     \item{\code{original}}{A \code{data.frame} of neighbors as generated
#'     by single-variable flips, before validation.}
#'     \item{\code{validated}}{A \code{data.frame} of neighbors after
#'     applying \code{validateModels()}, representing valid pharmacometric
#'     models.}
#'   }
#'
#' @examples
#' current_string <- c(no.cmpt = 2, eta.km = 0, eta.vc = 1,
#'                     eta.vp = 0, eta.vp2 = 0, eta.q = 1,
#'                     eta.q2 = 0, mm = 0, mcorr = 1, rv = 2)
#' neighbors <- generate_neighbors_df(current_string, search.space = "ivbase")
#' head(neighbors$original)   # raw neighbors (pre-validation)
#' head(neighbors$validated)  # validated neighbors (post-validation)
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

  # Store both original and validated neighbors
  orig_neighbors <- list()
  valid_neighbors <- list()

  # Generate neighbors by modifying one variable at a time
  for (variable in names(current_string)) {
    options <- options_list[[variable]]
    for (option in options) {
      if (option != current_string[variable]) {
        new_neighbor <- current_string
        new_neighbor[variable] <- option
        orig_neighbors <- append(orig_neighbors, list(new_neighbor))

        # validate this neighbor
        validated <- validateModels(
          string = as.numeric(new_neighbor),
          search.space = search.space,
          code.source = "TS"
        )
        valid_neighbors <- append(valid_neighbors, list(validated))
      }
    }
  }

  # Convert to data.frames
  orig_df  <- do.call(rbind, orig_neighbors)
  valid_df <- do.call(rbind, valid_neighbors)

  colnames(orig_df)  <- names(current_string)
  colnames(valid_df) <- names(current_string)

  # Deduplicate after validation (based on validated version)
  keep_idx <- !duplicated(valid_df)
  orig_df  <- orig_df[keep_idx, , drop = FALSE]
  valid_df <- valid_df[keep_idx, , drop = FALSE]

  # Candidate list sampling (optional)
  if (!is.null(candidate.size) &&
      nrow(valid_df) > candidate.size) {
    if (!is.null(seed))
      set.seed(seed)
    idx <- sample(1:nrow(valid_df), candidate.size)
    orig_df  <- orig_df[idx, , drop = FALSE]
    valid_df <- valid_df[idx, , drop = FALSE]
  }


  # Return both original and validated neighbors together
  return(list(
    original_neighbors  = as.data.frame(orig_df),
    validated_neighbors = as.data.frame(valid_df)
  ))
}



#' Detect the primary move between two model codes
#'
#' Compares a previous model code with a new one and identifies
#' the primary intended change. If an \code{original_neighbor}
#' is provided, this is used to determine the intended change,
#' ignoring any secondary modifications introduced by validation.
#'
#' @param prev_string A named numeric vector: the starting model code.
#' @param new_string A named numeric vector: the validated model code.
#' @param original_neighbor Optional named numeric vector: the original
#'   neighbor before validation. If provided, this is used to identify
#'   the primary change.
#'
#' @return A list with \code{element}, \code{from}, and \code{to}
#'   describing the primary change.
#'
#' @examples
#' prev <- c(no.cmpt = 2, eta.vc = 1)
#' orig <- c(no.cmpt = 3, eta.vc = 1)  # original neighbor
#' new  <- c(no.cmpt = 3, eta.vc = 0)  # validated neighbor (extra fix)
#' detect_move(prev, new, original_neighbor = orig)
#'
#' @export
detect_move <- function(prev_string, new_string, original_neighbor = NULL) {
  if (!is.null(original_neighbor)) {
    # Use original_neighbor to locate the intended flip
    diff_idx <- which(prev_string != original_neighbor)
  } else {
    # Fallback: detect based on new_string
    diff_idx <- which(prev_string != new_string)
  }

  if (length(diff_idx) == 0) {
    return(NULL)  # no change
  }
  if (length(diff_idx) > 1) {
    stop("Multiple changes detected; using the first difference as primary.")
    diff_idx <- diff_idx[1]
  }

  element <- names(prev_string)[diff_idx]
  from <- prev_string[diff_idx]
  to   <- new_string[diff_idx]

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
#' @param tabu.policy Character scalar. Tabu restriction type:
#'   \code{"attribute"} (default) or \code{"move"}.
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
is_move_tabu <- function(move, tabu_list, tabu.policy = c("attribute", "move")) {
  tabu.policy <- match.arg(tabu.policy)

  if (is.null(move) || is.null(tabu_list) || nrow(tabu_list) == 0) {
    return(FALSE)
  }

  if (tabu.policy == "move") {
    # Move-based tabu: forbid only the exact from→to transition
    return(any(
      tabu_list$element == move$element &
        tabu_list$from   == move$from &
        tabu_list$to     == move$to &
        tabu_list$tabu.iteration.left > 0
    ))
  } else if (tabu.policy == "attribute") {
    # Attribute-based tabu: forbid any move that sets the element to a tabu value
    return(any(
      tabu_list$element == move$element &
        tabu_list$to     == move$to &
        tabu_list$tabu.iteration.left > 0
    ))
  }
}

#' Apply 2-bit perturbation to escape local optimum
#'
#' This function randomly flips two parameters ("2-bit change") in the current
#' model string to generate a perturbed candidate. The candidate is then passed
#' through \code{validateModels()} to ensure pharmacological validity.
#'
#' The function returns both:
#' \itemize{
#'   \item \code{original_neighbor}: the raw 2-bit flip before validation
#'   \item \code{validated_neighbor}: the corrected version after validation
#' }
#' This allows downstream functions (e.g. \code{detect_move()}) to identify
#' which parameters were intentionally changed (primary moves), while still using
#' a valid model code for evaluation.
#'
#' @param prev_string A named numeric vector representing the current model.
#' @param search.space The search space ("ivbase" or "oralbase").
#' @param max.try Maximum number of attempts to generate a valid perturbed model.
#'
#' @return A \code{list} with two named numeric vectors:
#' \item{original_neighbor}{raw 2-bit flip (may be invalid)}
#' \item{validated_neighbor}{validated and usable model code}
#'
#' @examples
#' prev <- c(no.cmpt = 2, eta.km = 0, eta.vc = 1,
#'           eta.vp = 0, eta.vp2 = 0, eta.q = 1,
#'           eta.q2 = 0, mm = 0, mcorr = 1, rv = 2)
#' perturb <- perturb_2bit(prev, search.space = "ivbase")
#' perturb$original_neighbor   # original 2-bit flip
#' perturb$validated_neighbor  # validated model
#'
#' @export
perturb_2bit <-
  function(prev_string, search.space, max.try = 1000) {
    for (i in 1:max.try) {
      disturbed <- prev_string
      flip_idx <- sample(seq_along(prev_string), 2)  # flip 2 bits

      for (idx in flip_idx) {
        param <- names(prev_string)[idx]
        if (param %in% c("no.cmpt", "rv")) {
          all_options <- 1:3
        } else {
          all_options <- 0:1
        }

        current_val <- as.numeric(disturbed[idx])
        candidates <- setdiff(all_options, current_val)

        if (length(candidates) > 0) {
          disturbed[idx] <- sample(candidates, 1)
        }
      }

      # Save original 2-bit flip (before validation)
      disturbed_orig <- disturbed

      # Validate neighbor
      disturbed_val <- validateModels(string       = disturbed,
                                      search.space = search.space,
                                      code.source  = "TS")
      names(disturbed_val) <- names(prev_string)

      # Check that at least 2 positions changed after validation
      changed_positions <- sum(disturbed_val != prev_string)
      if (changed_positions >= 2) {
        return(list(
          original_neighbor  = disturbed_orig,
          validated_neighbor = disturbed_val
        ))
      }
    }
    warning("2-bit perturbation failed after ",
            max.try,
            " attempts, returning original string.")
    return(list(
      original_neighbor  = prev_string,
      validated_neighbor = prev_string
    ))
  }
