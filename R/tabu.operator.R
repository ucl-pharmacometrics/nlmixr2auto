#' Control Parameters for Tabu Search
#'
#' Creates a list of control settings for the \code{\link{tabu.operator}} function.
#'
#' @param tabu.duration Integer. Number of iterations a move remains tabu.
#' @param max.round Integer. Maximum number of search iterations.
#' @param start.point Optional initial model code vector. If NULL, defaults to a 2-compartment model.
#' @param aspiration Logical. Whether to apply the aspiration criterion.
#'   If TRUE, tabu moves are allowed if they yield a solution strictly better
#'   than the global best found so far.
#' @param candidate.size Optional integer. If not NULL, restricts neighborhood sear
#'   to a random subset of this size (candidate list strategy).
#' @param tabu.policy Character. Type of tabu restriction:
#'   \itemize{
#'     \item \code{"attribute"} — forbid revisiting a variable value (default).
#'     \item \code{"move"} — forbid only specific from–to transitions.
#'   }
#'
#' @return A list of Tabu Search hyperparameters.
#' @export
tabuControl <- function(tabu.duration = 3,
                        max.round = 20,
                        start.point = NULL,
                        aspiration = TRUE,
                        candidate.size = NULL,
                        tabu.policy="attribute") {
  list(
    tabu.duration = tabu.duration,
    max.round = max.round,
    start.point = start.point,
    aspiration = aspiration,
    candidate.size = candidate.size,
    tabu.policy=tabu.policy
  )
}

#' Tabu Search Operator for Pharmacometric Model Selection
#'
#' Performs Tabu Search to explore the pharmacometric model space and
#' identify the best-performing model. Supports both IV and Oral search spaces.
#'
#' @param dat Dataset (typically PK/PD data).
#' @param param_table Optional parameter table. If \code{NULL}, generated via \code{auto_param_table()}.
#' @param tabu.control A list of Tabu Search control parameters:
#'   \itemize{
#'     \item \code{tabu.duration} — Integer. Number of iterations a move remains tabu.
#'     \item \code{max.round} — Integer. Maximum number of search iterations.
#'     \item \code{start.point} — Optional initial model code vector.
#'     \item \code{aspiration} — Logical. If \code{TRUE}, allows aspiration criterion.
#'     \item \code{candidate.size} — Optional integer. Maximum number of neighbors
#'           randomly sampled from the full neighborhood (candidate list strategy).
#'   }
#' @param search.space Character. Search space type: \code{"ivbase"} or \code{"oralbase"}.
#' @param no.cores Integer. Number of CPU cores to use.
#' @param foldername Character. Folder name for temporary outputs.
#' @param filename Character. Base filename for outputs.
#' @param penalty.control A list of penalty settings, typically from \code{penaltyControl()}.
#' @param precomputed_results_file Optional cache file for model results.
#' @param seed.no Random seed for reproducibility.
#' @param ... Additional arguments passed to \code{mod.run()}.
#'
#' @return An object of class \code{tabuOperatorResult}, containing:
#'   \item{\code{Final Selected Code}}{Vector representation of the best model.}
#'   \item{\code{Final Selected Model Name}}{Selected best model (human-readable).}
#'   \item{\code{Model Run History}}{Data frame of all model evaluations with fitness values.}
#'   \item{\code{Search History}}{List with iteration-level history:
#'         \code{starting.points.history}, \code{local.best.history},
#'         \code{tabu.elements.history}, \code{neighbors.history}.}
#'
#' @details
#' This function implements a customized Tabu Search framework for pharmacometric
#' model structure optimization. Key design aspects:
#'
#' - **Solution representation**: Models are encoded as bit vectors
#'   (e.g., number of compartments, inclusion of random effects, residual error structure).
#'
#' - **Initial solution**: Default is a 2-compartment base model, or user-specified
#'   via \code{tabu.control$start.point}.
#'
#' - **Neighborhood definition**: For each iteration, all one-bit flips are generated
#'   (\emph{original neighbors}), then passed through \code{validateModels()} to ensure
#'   only valid pharmacometric models are retained (\emph{validated neighbors}).
#'
#' - **Move definition**: A move is defined based on the \emph{original neighbor}
#'   (primary 1-bit flip), even if the validated neighbor differs after correction.
#'
#' - **Tabu list**: Stores recent moves with remaining tabu tenure (\code{tabu.duration}).
#'
#' - **Aspiration criterion**: If enabled, a tabu move is allowed if it produces
#'   a model better than the historical best.
#'
#' - **Perturbation**: If no improving neighbor is found, a random 1-bit perturbation
#'   is applied. Both original and validated perturbed neighbors are tracked.
#'
#' - **Objective function**: Model fit quality (e.g., AIC/OFV), computed by \code{mod.run()},
#'   with penalties applied via \code{penalty.control}.
#'
#' - **Termination**: Fixed number of iterations (\code{max.round}).
#'
#' - **History tracking**: For each iteration, the algorithm stores:
#'   starting points, validated neighbors, local best models, and tabu elements.
#'
#' @seealso \code{\link{generate_neighbors_df}}, \code{\link{detect_move}},
#'          \code{\link{perturb_1bit}}, \code{\link{penaltyControl}},
#'          \code{\link{CodetoMod}}
#'
#' @export

tabu.operator <- function(dat,
                          param_table = NULL,
                          tabu.control = tabuControl(),
                          search.space = "ivbase",
                          no.cores = rxode2::getRxThreads(),
                          foldername = "test",
                          filename = "test",
                          penalty.control = penaltyControl(),
                          precomputed_results_file = NULL,
                          seed.no = 1234,
                          ...) {
  current.date <- Sys.Date()
  set.seed(seed.no)
  # tabu control
  tabu.duration <- tabu.control$tabu.duration
  max.round <- tabu.control$max.round
  start.point <- tabu.control$start.point
  aspiration <- tabu.control$aspiration
  candidate.size <- tabu.control$candidate.size
  tabu.policy <- tabu.control$tabu.policy

  # Create temporary output directory
  outputdir <- paste0("TS_",
                      current.date,
                      "-",
                      foldername,
                      "_",
                      digest::digest(dat),
                      "_temp")
  if (!dir.exists(outputdir)) {
    dir.create(outputdir, showWarnings = FALSE, recursive = TRUE)
  }
  setwd(file.path(getwd(), outputdir))
  storage.path <- getwd()

  # --- Initialize parameter table ---
  param_table <- auto_param_table(
    dat = dat,
    param_table = param_table,
    nlmixr2autoinits = TRUE,
    foldername = foldername
  )

  # --- Define bit names for search space ---
  if (search.space == "ivbase") {
    bit.names <- c(
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
  } else if (search.space == "oralbase") {
    bit.names <- c(
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
    )
  } else {
    stop("Unknown search.space type: must be 'ivbase' or 'oralbase'")
  }

  # --- Initialize starting point ---
  if (is.null(tabu.control$start.point)) {
    string_vec <- base_model(search.space)
    string_vec["no.cmpt"] <- 2 # enforce 2-compartment default
  } else {
    string_vec <- tabu.control$start.point
  }

  # --- Initialize histories ---
  local.best.history <- list()
  starting.points.history <- list()
  tabu.elements.history <- list()
  neighbors.history <- list()
  tabu.elements.all <- NULL
  prev_string <- NULL
  globalbest <- NULL

  pb <-
    progress::progress_bar$new(
      format = " Tabu Search [:bar] :percent (iteration :current/:total)\n",
      total = max.round,
      clear = FALSE,
      width = 60
    )

    for (r in 1:max.round) {
      # Step 1: define current starting point
      if (r == 1) {
        start_string <- string_vec
      } else {
        start_string <- current_string
      }
      starting.points.history[[r]] <- start_string

      # Step 2: generate neighbors (original + validated)
      neighbors_list <-
        generate_neighbors_df(start_string, search.space = search.space)
      neighbors_orig <-
        neighbors_list$original_neighbors   # pre-validation
      neighbors_val  <-
        neighbors_list$validated_neighbors  # post-validation (legal models)

      # Deduplicate validated neighbors
      neighbors_val <- distinct(neighbors_val)[, bit.names]

      # Step 3a: Filter tabu neighbors (without aspiration logic)
      neighbors_eval <- list()
      aspiration_candidates <- list()
      move_records <- list()  # store primary moves

      if (!is.null(tabu.elements.all) &&
          nrow(tabu.elements.all) > 0) {
        for (row in 1:nrow(neighbors_val)) {
          # detect primary move using original neighbor
          move <- detect_move(start_string,
                              neighbors_val[row, ],
                              original_neighbor = neighbors_orig[row, ])

          if (is_move_tabu(move = move, tabu_list = tabu.elements.all,tabu.policy = tabu.policy)) {
            # tabu move
            if (tabu.control$aspiration) {
              aspiration_candidates[[length(aspiration_candidates) + 1]] <-
                neighbors_val[row, ]
            }
          } else {
            # non-tabu -> keep
            neighbors_eval[[length(neighbors_eval) + 1]] <-
              neighbors_val[row, ]
            move_records[[length(move_records) + 1]] <- move$element
          }
        }
      } else {
        # if no tabu elements, all neighbors are valid
        neighbors_eval <-
          split(neighbors_val, seq_len(nrow(neighbors_val)))
        # detect moves for all neighbors
        move_records <-
          lapply(seq_len(nrow(neighbors_val)), function(row) {
            move <- detect_move(start_string,
                                neighbors_val[row, ],
                                original_neighbor = neighbors_orig[row, ])
            move$element
          })
      }

      # Convert lists back to data frames
      if (length(neighbors_eval) > 0) {
        neighbors_eval <- do.call(rbind, neighbors_eval)
        neighbors_eval$move.element <- unlist(move_records)
      } else {
        neighbors_eval <- data.frame()
      }

      if (length(aspiration_candidates) > 0) {
        aspiration_candidates <- do.call(rbind, aspiration_candidates)
      } else {
        aspiration_candidates <- data.frame()
      }

      # save (store validated neighbors in history)
      neighbors.history[[r]] <- neighbors_val

      # --- Step 4: Evaluate neighbors (fitness calculation) ---
      if (nrow(neighbors_eval) > 0) {
        neighbors_eval$fitness <- vapply(seq_len(nrow(neighbors_eval)),
                                         function(k) {
                                           string_vec <- as.numeric(neighbors_eval[k, bit.names])
                                           result <- try(mod.run(
                                             r = r,
                                             dat = dat,
                                             search.space = search.space,
                                             string = string_vec,
                                             param_table = param_table,
                                             penalty.control = penalty.control,
                                             precomputed_results_file = precomputed_results_file,
                                             filename = filename
                                           ),
                                           silent = TRUE)
                                           if (is.numeric(result) &&
                                               length(result) == 1) {
                                             result
                                           } else {
                                             NA_real_
                                           }
                                         },
                                         numeric(1))

        # --- Aspiration criterion check ---
        if (tabu.control$aspiration) {
          best.fitness <- min(Store.all$fitness, na.rm = TRUE)
          aspiration_candidates <-
            neighbors_eval[neighbors_eval$fitness < best.fitness, ]
        } else {
          aspiration_candidates <- data.frame()
        }
      } else {
        neighbors_eval <- data.frame()
        aspiration_candidates <- data.frame()
      }

      # Step 5: update local best
      if (nrow(neighbors_eval) > 0) {
        localbest <-
          neighbors_eval[which.min(neighbors_eval$fitness), , drop = FALSE]
      } else {
        localbest <- start_string  # fallback, in case no neighbors
      }
      local.best.history[[r]] <- localbest

      if (r == 1) {
        prev_string <- string_vec      # starting point
      } else {
        prev_string <- start_string
      }
      current_string <- localbest[1, bit.names]

      # --- Step X: Update current solution with neighbor or perturbation ---
      has_been_start <- any(vapply(starting.points.history, function(hist) {
        all(hist == current_string)
      }, logical(1)))

      if (has_been_start) {
        message(
          "Iteration ", r,
          ": candidate already used as a starting point. Applying 2-bit perturbation to avoid cycling."
        )
        perturb <- perturb_2bit(prev_string, search.space)
        current_string <- perturb$validated_neighbor

        tabu.elements <- data.frame(
          tabu.num = r,
          element  = "perturbation", # not in the tabulist
          from     = NA,
          to       = NA,
          tabu.iteration.left = 0,
          stringsAsFactors = FALSE
        )
      } else {
        # Case 2: Normal neighbor move update tabu as usual
        idx <- match(paste0(as.numeric(current_string), collapse = "_"),
                     apply(neighbors_val[, bit.names], 1, function(x)
                       paste0(as.numeric(x), collapse = "_")))

        if (!is.na(idx)) {
          # Found the matching original neighbor → record the true move
          move <- detect_move(start_string,
                              new_string        = current_string,
                              original_neighbor = neighbors_orig[idx,])
        } else {
          # Fallback: in rare cases where match fails,
          move <- detect_move(start_string,
                              new_string        = current_string,
                              original_neighbor = current_string
          )
        }

        if (tabu.policy == "move") {
          # Move-based tabu:
          # Store both forward and reverse moves (e.g., 2 to 3 and 3 to 2)
          tabu.elements <- rbind(
            data.frame(
              tabu.num = r,
              element  = move$element,
              from     = unname(move$from),
              to       = unname(move$to),
              tabu.iteration.left = tabu.control$tabu.duration,
              stringsAsFactors = FALSE
            ),
            data.frame(
              tabu.num = r,
              element  = move$element,
              from     = unname(move$to),   # reverse move
              to       = unname(move$from), # reverse move
              tabu.iteration.left = tabu.control$tabu.duration,
              stringsAsFactors = FALSE
            )
          )
        } else if (tabu.policy == "attribute") {
          # Attribute-based tabu:
          # Store only the target value (e.g., "element = no.cmpt, to = 3")
          # This forbids any move that sets the element to this value.
          tabu.elements <- data.frame(
            tabu.num = r,
            element  = move$element,
            from     = unname(move$from),
            to       = unname(move$to),
            tabu.iteration.left = tabu.control$tabu.duration,
            stringsAsFactors = FALSE
          )
        }
      }

      if (!is.null(tabu.elements.all)) {
        tabu.elements.all$tabu.iteration.left <-
          tabu.elements.all$tabu.iteration.left - 1
      }

      tabu.elements.all <- rbind(tabu.elements.all, tabu.elements)
      tabu.elements.all <-
        tabu.elements.all[tabu.elements.all$tabu.iteration.left > 0, ]

      rownames(tabu.elements.all) <- NULL
      tabu.elements.history[[r]] <- tabu.elements.all

      pb$tick()
    } # end loop

  # ----------------------------
  # Final output (Tabu Search)
  # ----------------------------
  localbestf <- Store.all[Store.all$fitness == min(Store.all$fitness), ][1, ]
  best_model_code <- as.numeric(localbestf[, bit.names])
  names(best_model_code) <- bit.names
  best_model_name <- CodetoMod(sel.best.code = best_model_code,
                               search.space  = search.space)

  out <- new.env(parent = emptyenv())
  class(out) <- "tabuOperatorResult"
  out[["Final Selected Code"]] <- best_model_code
  out[["Final Selected Model Name"]] <- best_model_name
  out[["Model Run History"]] <- as.data.frame(Store.all, stringsAsFactors = FALSE)
  out[["Search History"]] <- list(
    starting.points.history = starting.points.history,
    local.best.history      = local.best.history,
    tabu.elements.history   = tabu.elements.history,
    neighbors.history       = neighbors.history
  )

  on.exit({
    rm(modi, r, Store.all, precomputed_cache_loaded, envir = .GlobalEnv)
  }, add = TRUE)

  return(out)

}


#' Print Method for Tabu Search Results
#'
#' @description
#' This function defines the print method for objects of class
#' \code{tabuOperatorResult}, which are returned by the
#' \code{\link{tabu.operator}} function. It prints the final selected
#' model code and the corresponding human-readable model name obtained
#' by the Tabu Search algorithm.
#'
#' @param x An object of class \code{tabuOperatorResult}, typically the
#'   output from \code{\link{tabu.operator}}.
#' @param ... Additional arguments (currently not used).
#'
#' @details
#' The function prints:
#' - The **final selected model code**, as a named numeric vector.
#' - The **final selected model name**, generated via \code{CodetoMod()}.
#'
#' Colored console output is provided using the \pkg{crayon} package,
#' to visually distinguish results from other optimization operators
#' (e.g., ACO).
#'
#' @return
#' The function is called for its side effects (printing). It returns
#' the input object \code{x} invisibly.
#'
#' @seealso
#' \code{\link{tabu.operator}}, \code{\link{CodetoMod}}
#'
#' @examples
#' \dontrun{
#' result <- tabu.operator(dat = mydata, search.space = "ivbase")
#' print(result)   # automatically calls print.tabuOperatorResult
#' }
#'
#' @export
print.tabuOperatorResult <- function(x, ...) {
  # Print final selected model code
  cat(crayon::green$bold("\n=== Final Selected Model Code (Tabu Search) ===\n"))
  print(x$`Final Selected Code`)

  # Print final selected model name
  cat(crayon::green$bold("\n=== Final Selected Model Name (Tabu Search) ===\n"))
  cat(x$`Final Selected Model Name`, "\n")

  invisible(x)
}



