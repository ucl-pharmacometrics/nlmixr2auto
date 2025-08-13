#' Create a base model code for stepwise model search
#'
#' Generates a named numeric vector representing the base model code
#' for a given search space (\code{"ivbase"} or \code{"oralbase"}).
#' The returned vector contains model specification fields such as
#' number of compartments, IIV flags, Michaelis–Menten term, correlation flag,
#' and residual error model.
#'
#' Users can optionally supply a custom base model code via \code{custom_base}.
#' The function will validate its type and length according to the chosen
#' \code{search.space} and return it with proper element names.
#'
#' @param search.space Character scalar: either \code{"ivbase"} (9 elements)
#'   or \code{"oralbase"} (10 elements). Default is \code{"ivbase"}.
#' @param custom_base Optional numeric vector representing a user-specified
#'   base model code. Must have length 9 for \code{"ivbase"} or 10 for
#'   \code{"oralbase"}.
#'
#' @return A named integer vector representing the model code, with elements:
#'   \itemize{
#'     \item \code{no.cmpt}   — Number of compartments
#'     \item \code{eta.km}    — IIV flag for Km
#'     \item \code{eta.vc}    — IIV flag for Vc
#'     \item \code{eta.vp}    — IIV flag for Vp
#'     \item \code{eta.vp2}   — IIV flag for Vp2
#'     \item \code{eta.q}     — IIV flag for Q
#'     \item \code{eta.q2}    — IIV flag for Q2
#'     \item \code{mm}        — Michaelis–Menten term flag
#'     \item \code{mcorr}     — Correlation flag between ETAs
#'     \item \code{eps.model} — Residual error model code
#'   }
#'   For \code{"ivbase"}, the last element (\code{eps.model}) is in position 9.
#'
#' @examples
#' base_model("ivbase")
#' base_model("oralbase")
#' base_model("ivbase", custom_base = c(2,1,0,0,0,0,0,0,3))
#'
#' @export
base_model <- function(search.space = "ivbase",
                       custom_base = NULL) {
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

  if (!is.null(custom_base)) {
    if (!is.numeric(custom_base)) {
      stop("custom_base must be a numeric vector.")
    }
    if (search.space == "ivbase" && length(custom_base) != 9L) {
      stop("For search.space = 'ivbase', model code must have length 9.")
    }
    if (search.space == "oralbase" && length(custom_base) != 10L) {
      stop("For search.space = 'oralbase', model code must have length 10.")
    }
    return(setNames(
      as.integer(custom_base),
      if (search.space == "ivbase")
        element_names[1:9]
      else
        element_names
    ))
  }

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






#' Screen number of compartments from the current/baseline model
#'
#' Clones the current model code (from \code{state$best_code}, or from
#' \code{base_model(search.space, custom_base = ...)} if \code{state$best_code}
#' is \code{NULL}) and generates three candidates by setting
#' \code{no.cmpt} to 1, 2, and 3 while keeping all other positions unchanged.
#' Each candidate is evaluated with \code{mod.run()} and the best (lowest
#' \code{Fitness}) is returned alongside a results table for logging.
#'
#' @param dat A pharmacokinetic dataset passed through to \code{mod.run()}.
#' @param state A list-like state that may contain:
#'   \itemize{
#'     \item \code{best_code}: named integer vector for the current model code (may be \code{NULL})
#'     \item \code{modi}: optional value forwarded to \code{mod.run()} (if absent, defaults to 1)
#'   }
#' @param search.space Character scalar, e.g. \code{"ivbase"} or \code{"oralbase"}.
#' @param param_table Object passed through to \code{mod.run()}.
#' @param ... Additional arguments forwarded to \code{mod.run()}. Must include
#'   \code{custom_base} (a numeric vector) used by \code{base_model()} to
#'   initialize the starting code when \code{state$best_code} is \code{NULL}.
#'
#' @details
#' Candidates are formed by copying the starting code and changing only the
#' \code{no.cmpt} element to \{1, 2, 3\}. Named codes are converted to unnamed
#' numeric vectors (\code{unname()}) when passed to \code{mod.run()}.
#' If \code{state$modi} is \code{NULL}, a default of \code{1} is used.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{results_table}: \code{data.frame} with columns
#'         \code{Step}, \code{"Model name"}, \code{"Model code"}, \code{Fitness}
#'   \item \code{best_code}: named integer vector of the best candidate’s model code
#'   \item \code{best_row}: one-row \code{data.frame} (the best candidate)
#' }
#'
#' @examples
#' # Initialize state (if needed)

#' # Run the compartment screening (custom_base is required and forwarded via ...):
#'res <- step_compartments(
#'  dat = dat,
#'  state = state,
#'  search.space = "ivbase",
#'  param_table = param_table,
#'  saem.control = saemControl(nBurn = 10,nEm = 10),
#'  table.control = tableControl(cwres = T)
#')
#'
#' # Log both all candidates and the best candidate, and update state$best_code:
#' state <- modlog_state(state, results_table = res$results_table, step_name = "no. of compartments")
#'
#' @export

step_compartments <-
  function(dat, state, search.space, param_table,penalty.control=NULL,precomputed_results_file=NULL,filename=NULL,...) {
    dots <- list(...)
    # custom_base <- dots$custom_base
    if (!is.null(dots$custom_base)) {
      custom_base <- dots$custom_base
    } else {
      custom_base <- NULL
    }

    if (!is.null(state$best_code)) {
      current_code <- state$best_code
    } else {
      current_code <- base_model(search.space=search.space, custom_base = custom_base)
    }

    candidate_codes <- lapply(1:3, function(k) {
      code <- current_code
      code["no.cmpt"] <- k
      code
    })

    fits <- vapply(candidate_codes, function(code_vec) {
      mod.run(
        modi         = modi,
        string       = unname(code_vec),
        dat          = dat,
        search.space = search.space,
        param_table  = param_table,
        precomputed_results_file =   precomputed_results_file,
        filename=filename,
        ...
      )
    }, numeric(1))

    # 4) results table
    model_names <- vapply(candidate_codes,
                          function(code)
                            CodetoMod(search.space, unname(code)),
                          character(1))
    model_codes_chr <- vapply(candidate_codes,
                              function(code)
                                paste(unname(code), collapse = ","),
                              character(1))

    results_table <- data.frame(
      Step          = "No. of compartments",
      "Penalty terms" = paste(penalty.control$penalty.terms, collapse = ", "),
      "Model name"  = model_names,
      "Model code"  = model_codes_chr,
      Fitness       = fits,
      stringsAsFactors = FALSE
    )

    # 5) pick best
    best_idx  <- which.min(results_table$Fitness)
    best_row  <- results_table[best_idx, , drop = FALSE]
    best_code <- candidate_codes[[best_idx]]

    r<<-r+1
    list(results_table = results_table,
         best_code     = best_code,
         best_row      = best_row)
  }


#' Screen elimination type (Michaelis–Menten on/off)
#'
#' This step compares a standard linear elimination model (`mm = 0`)
#' versus a Michaelis–Menten elimination model (`mm = 1`) while keeping
#' all other positions of the current model code unchanged.
#'
#' If `mm = 0`, any IIV term for Km (`eta.km`) will be automatically set to 0.
#'
#' @param dat Dataset passed to \code{mod.run()}.
#' @param state List-like with optional fields:
#'   \itemize{
#'     \item \code{best_code}: named integer vector (starting code). If NULL, \code{base_model()} is used.
#'     \item \code{modi}: forwarded to \code{mod.run()} (default \code{1} if NULL).
#'   }
#' @param search.space Character scalar: either \code{"ivbase"} or \code{"oralbase"}.
#' @param param_table Parameter table passed to \code{mod.run()}.
#' @param ... Additional arguments forwarded to \code{mod.run()}. May include
#'   \code{custom_base} used by \code{base_model()}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{results_table}: \code{data.frame} with columns Step, Model name, Model code, Fitness
#'   \item \code{best_code}: named integer vector of the best candidate’s model code
#'   \item \code{best_row}: one-row \code{data.frame} (the best candidate)
#' }
#'
#' @examples
#' \dontrun{
#' # Example dataset and parameter table would normally come from your workflow
#' dat <- your_data_frame
#' param_table <- your_param_table
#'
#' # Initialize state with a base model
#' state <- list(best_code = base_model("ivbase"), modi = 1)
#'
#' # Run elimination type screening
#' res_mm <- step_elimination(
#'   dat = dat,
#'   state = state,
#'   search.space = "ivbase",
#'   param_table = param_table
#' )
#'
#' # View results
#' res_mm$results_table
#' res_mm$best_code
#' }
#'
#' @export
step_elimination <-
  function(dat, state, search.space, param_table,penalty.control=NULL,precomputed_results_file=NULL,filename=NULL, ...) {
    dots <- list(...)
    if (!is.null(dots$custom_base)) {
      custom_base <- dots$custom_base
    } else {
      custom_base <- NULL
    }
    # starting code
    if (!is.null(state$best_code)) {
      current_code <- state$best_code
    } else {
      current_code <-
        base_model(search.space = search.space, custom_base = custom_base)
    }

    # candidates: mm = 0 and mm = 1
    candidate_codes <- lapply(c(0L, 1L), function(m) {
      code <- current_code
      code["mm"] <- m
      # Optional: if no MM, Km IIV should be off
      if (m == 0L && "eta.km" %in% names(code))
        code["eta.km"] <- 0L
      code
    })

    fits <- vapply(candidate_codes, function(code_vec) {
      mod.run(
        modi         = modi,
        string       = unname(code_vec),
        dat          = dat,
        search.space = search.space,
        param_table  = param_table,
        precomputed_results_file=   precomputed_results_file,
        filename = filename,
        ...
      )
    }, numeric(1))

    to_name <- function(code) {
      if (exists("CodetoMod", mode = "function")) {
        CodetoMod(search.space, unname(code))
      } else {
        paste0("Model_", paste(unname(code), collapse = "_"))
      }
    }

    model_names <- vapply(candidate_codes, to_name, character(1))
    model_codes_chr <-
      vapply(candidate_codes, function(code)
        paste(unname(code), collapse = ","), character(1))

    results_table <- data.frame(
       Step        = "Elimination type",
       "Penalty terms" = paste(penalty.control$penalty.terms, collapse = ", "),
      "Model name" = model_names,
      "Model code" = model_codes_chr,
      Fitness      = fits,
      stringsAsFactors = FALSE
    )

    best_idx <- which.min(results_table$Fitness)

    r<<-r+1

    list(
      results_table = results_table,
      best_code     = candidate_codes[[best_idx]],
      best_row      = results_table[best_idx, , drop = FALSE]
    )
  }



#' Stepwise inclusion of inter-individual variability (IIV)
#'
#' Two-stage procedure:
#' 1) Special IIV:
#'    - If \code{mm = 1} and \code{eta.km = 0}, test \code{eta.km = 1}.
#'    - If \code{search.space = "oralbase"} and \code{eta.ka = 0}, test \code{eta.ka = 1}.
#' 2) Structural IIV (forward selection):
#'    - 1-compartment: \code{eta.vc}
#'    - 2-compartment: \code{eta.vc, eta.vp, eta.q}
#'    - 3-compartment: \code{eta.vc, eta.vp, eta.q, eta.vp2, eta.q2}
#'    At each iteration, turn exactly one available IIV from 0 to 1, evaluate all candidates,
#'    accept the single best improvement, and stop when no improvement occurs.
#'
#' @param dat Data frame passed to \code{mod.run()}.
#' @param state List-like object with:
#'   - \code{best_code}: starting named integer/numeric vector; if \code{NULL}, \code{base_model()} is used.
#'   - \code{modi}: optional; forwarded to \code{mod.run()} (default \code{1L}).
#' @param search.space Character scalar: \code{"ivbase"} or \code{"oralbase"}.
#' @param param_table Parameter table passed to \code{mod.run()}.
#' @param penalty.control Optional penalty control object passed to \code{mod.run()}.
#' @param ... Additional arguments forwarded to \code{mod.run()}. May include
#'   \code{custom_base} for \code{base_model()} when \code{state$best_code} is \code{NULL}.
#'
#' @return A list with:
#'   - \code{results_table}: \code{data.frame} with columns Step, Model name, Model code, Fitness.
#'   - \code{best_code}: final best model code (named integer vector).
#'   - \code{best_row}: one-row \code{data.frame} for the best entry.
#'
#' @examples
#' \dontrun{
#'
#' dat <- pheno_sd
#' param_table <-   param_table<-auto_param_table(dat = dat,
#'                                 param_table = param_table,
#'                                 nlmixr2autoinits = T)
#'results_iiv <- step_iiv(
#'  dat = dat,
#'  state = list(),
#'  search.space = "ivbase",
#'  param_table = param_table,
#'  penalty.control =   penalty.control
#')
#'
#' results_iiv$results_table
#' results_iiv$best_code
#'results_iiv$best_row
#' }
#'
#' @export
step_iiv <- function(dat,
                     state=list(),
                     search.space="ivbase",
                     param_table=NULL,
                     penalty.control = NULL,
                     precomputed_results_file=NULL,
                     filename=NULL,
                     ...) {

  dots <- list(...)
  if (!is.null(dots$custom_base)) {
    custom_base <- dots$custom_base
  } else {
    custom_base <- NULL
  }

  # Starting code
  if (!is.null(state$best_code)) {
    current_code <- state$best_code
  } else {
    current_code <- base_model(search.space = search.space, custom_base = custom_base)
  }

  all_results <- NULL

  # ---------------- Special IIV: eta.km (only if mm = 1) ----------------
  if ("mm" %in% names(current_code) && current_code["mm"] == 1L &&
      "eta.km" %in% names(current_code) && current_code["eta.km"] == 0L) {

    candidate_codes <- list(
      current_code,
      {tmp <- current_code; tmp["eta.km"] <- 1L; tmp}
    )

    fits <- vapply(candidate_codes, function(code_vec) {
      mod.run(
        modi             = modi,
        string           = unname(code_vec),
        dat              = dat,
        search.space     = search.space,
        param_table      = param_table,
        penalty.control  = penalty.control,
        precomputed_results_file =   precomputed_results_file,
        filename=filename,
        ...
      )
    }, numeric(1))

    model_names <- vapply(candidate_codes, function(code) {
      if (exists("CodetoMod", mode = "function")) {
        CodetoMod(search.space = search.space, sel.best.code = unname(code))
      } else {
        paste0("Model_", paste(unname(code), collapse = "_"))
      }
    }, character(1))

    model_codes_chr <- vapply(candidate_codes, function(code) {
      paste(unname(code), collapse = ",")
    }, character(1))

    results_table <- data.frame(
      Step         = "IIV on Km",
      "Penalty terms" = paste(penalty.control$penalty.terms, collapse = ", "),
      "Model name" = model_names,
      "Model code" = model_codes_chr,
      Fitness      = fits,
      stringsAsFactors = FALSE
    )
    all_results <- rbind(all_results, results_table)

    best_idx <- which.min(results_table$Fitness)
    current_code <- candidate_codes[[best_idx]]
  }

  # ---------------- Special IIV: eta.ka (only for oralbase) ----------------
  if (identical(search.space, "oralbase") &&
      "eta.ka" %in% names(current_code) && current_code["eta.ka"] == 0L) {

    candidate_codes <- list(
      current_code,
      {tmp <- current_code; tmp["eta.ka"] <- 1L; tmp}
    )

    fits <- vapply(candidate_codes, function(code_vec) {
      mod.run(
        modi             = modi,
        string           = unname(code_vec),
        dat              = dat,
        search.space     = search.space,
        param_table      = param_table,
        penalty.control  = penalty.control,
        precomputed_results_file=   precomputed_results_file,
        filename=filename,
        ...
      )
    }, numeric(1))

    model_names <- vapply(candidate_codes, function(code) {
      if (exists("CodetoMod", mode = "function")) {
        CodetoMod(search.space = search.space, sel.best.code = unname(code))
      } else {
        paste0("Model_", paste(unname(code), collapse = "_"))
      }
    }, character(1))

    model_codes_chr <- vapply(candidate_codes, function(code) {
      paste(unname(code), collapse = ",")
    }, character(1))

    results_table <- data.frame(
      Step         = "IIV on Ka",
      "Penalty terms" = paste(penalty.control$penalty.terms, collapse = ", "),
      "Model name" = model_names,
      "Model code" = model_codes_chr,
      Fitness      = fits,
      stringsAsFactors = FALSE
    )
    all_results <- rbind(all_results, results_table)

    best_idx <- which.min(results_table$Fitness)
    current_code <- candidate_codes[[best_idx]]
  }

  # ---------------- Structural IIV (forward selection) ----------------
  struct_iiv <- c("eta.vc")

  no.cmpt <- suppressWarnings(as.integer(current_code[["no.cmpt"]]))
  if (no.cmpt >= 2L) struct_iiv <- c(struct_iiv, "eta.vp", "eta.q")
  if (no.cmpt >= 3L) struct_iiv <- c(struct_iiv, "eta.vp2", "eta.q2")
  # struct_iiv <- struct_iiv[struct_iiv %in% names(current_code)]

  keep_going <- TRUE
  while (keep_going) {
    available <- struct_iiv[current_code[struct_iiv] == 0L]
    if (length(available) == 0L) break
    # Baseline + "add one IIV" candidates
    candidate_codes <- vector("list", length(available) + 1L)
    candidate_codes[[1L]] <- current_code
    for (i in seq_along(available)) {
      nm <- available[i]
      tmp <- current_code
      tmp[nm] <- 1L
      candidate_codes[[i + 1L]] <- tmp
    }

    fits <- vapply(candidate_codes, function(code_vec) {
      mod.run(
        modi             = modi,
        string           = unname(code_vec),
        dat              = dat,
        search.space     = search.space,
        param_table      = param_table,
        penalty.control  = penalty.control,
        precomputed_results_file=   precomputed_results_file,
        filename=filename,
        ...
      )
    }, numeric(1))

    model_names <- vapply(candidate_codes, function(code) {
      if (exists("CodetoMod", mode = "function")) {
        CodetoMod(search.space = search.space, sel.best.code = unname(code))
      } else {
        paste0("Model_", paste(unname(code), collapse = "_"))
      }
    }, character(1))

    model_codes_chr <- vapply(candidate_codes, function(code) {
      paste(unname(code), collapse = ",")
    }, character(1))

    results_table <- data.frame(
      Step         = "IIV (forward)",
      "Penalty terms" = paste(penalty.control$penalty.terms, collapse = ", "),
      "Model name" = model_names,
      "Model code" = model_codes_chr,
      Fitness      = fits,
      stringsAsFactors = FALSE
    )
    all_results <- rbind(all_results, results_table)

    best_idx <- which.min(results_table$Fitness)

    # If baseline remains best (index 1), stop; else accept improvement and continue
    if (best_idx == 1L) {
      keep_going <- FALSE
    } else {
      current_code <- candidate_codes[[best_idx]]
    }
  }

  # ---------------- Finalize outputs ----------------
  if (!is.null(all_results)) {
    best_row <- all_results[which.min(all_results$Fitness), , drop = FALSE]
  } else {

    base_fit <- mod.run(
      modi             = modi,
      string           = unname(current_code),
      dat              = dat,
      search.space     = search.space,
      param_table      = param_table,
      penalty.control  = penalty.control,
      precomputed_results_file=   precomputed_results_file,
      filename=filename,
      ...
    )
    best_row <- data.frame(
      Step = "IIV (none)",
      "Model name" = if (exists("CodetoMod", mode = "function")) {
        CodetoMod(search.space = search.space, sel.best.code = unname(current_code))
      } else {
        paste0("Model_", paste(unname(current_code), collapse = "_"))
      },
      "Model code" = paste(unname(current_code), collapse = ","),
      Fitness = base_fit,
      stringsAsFactors = FALSE
    )
    all_results <- best_row
  }
  r<<-r+1
  list(
    results_table = all_results,
    best_code     = current_code,
    best_row      = best_row
  )
}

#' Screen ETA correlation (on vs off)
#'
#' Tests whether enabling ETA correlation (\code{mcorr = 1}) improves the fitness
#' over keeping it disabled (\code{mcorr = 0}), while keeping all other fields fixed.
#'
#' @param dat Data frame passed to \code{mod.run()}.
#' @param state List-like with:
#'   - \code{best_code}: starting named integer/numeric vector; if \code{NULL}, \code{base_model()} is used.
#'   - \code{modi}: optional; forwarded to \code{mod.run()} (default \code{1L}).
#' @param search.space Character scalar: \code{"ivbase"} or \code{"oralbase"}.
#' @param param_table Parameter table passed to \code{mod.run()}.
#' @param penalty.control Optional penalty control object forwarded to \code{mod.run()}.
#' @param ... Additional arguments forwarded to \code{mod.run()}. May include
#'   \code{custom_base} for \code{base_model()} when \code{state$best_code} is \code{NULL}.
#'
#' @return A list with:
#'   - \code{results_table}: \code{data.frame} with columns Step, Model name, Model code, Fitness.
#'   - \code{best_code}: named integer vector of the best candidate’s model code.
#'   - \code{best_row}: one-row \code{data.frame} (the best candidate).
#'
#' @examples
#' \dontrun{
#' dat <- your_data
#' param_table <- your_param_table
#' state <- list(best_code = base_model("ivbase"), modi = 1L)
#'
#'result_corr <- step_correlation(
#'  dat = dat,
#'  search.space = "ivbase",
#'  param_table = param_table,
#'  penalty.control = penaltyControl()
#')
#'result_corr$results_table
#'result_corr$best_code
#'result_corr$best_row
#'}
#'
#' @export
step_correlation <- function(dat,
                             state=list(),
                             search.space="ivbase",
                             param_table=NULL,
                             penalty.control = NULL,
                             precomputed_results_file=NULL,
                             filename=NULL,
                             ...) {
  dots <- list(...)
  if (!is.null(dots$custom_base)) {
    custom_base <- dots$custom_base
  } else {
    custom_base <- NULL
  }

  # starting code
  if (!is.null(state$best_code)) {
    current_code <- state$best_code
  } else {
    current_code <-
      base_model(search.space = search.space, custom_base = custom_base)
  }

  corrcode <- suppressWarnings(as.integer(current_code["mcorr"]))
  alt_val <- if (corrcode == 0L)
    1L
  else
    0L

  candidate_codes <- list(current_code,
                          {
                            tmp <- current_code
                            tmp["mcorr"] <- alt_val
                            tmp
                          })

  fits <- vapply(candidate_codes, function(code_vec) {
    mod.run(
      modi             = modi,
      string           = unname(code_vec),
      dat              = dat,
      search.space     = search.space,
      param_table      = param_table,
      penalty.control  = penalty.control,
      precomputed_results_file= precomputed_results_file,
      filename=filename,
      ...
    )
  }, numeric(1))

  model_names <- vapply(candidate_codes, function(code) {
    if (exists("CodetoMod", mode = "function")) {
      CodetoMod(search.space = search.space,
                sel.best.code = unname(code))
    } else {
      paste0("Model_", paste(unname(code), collapse = "_"))
    }
  }, character(1))

  model_codes_chr <- vapply(candidate_codes, function(code) {
    paste(unname(code), collapse = ",")
  }, character(1))

  results_table <- data.frame(
    Step         = "ETA correlation",
    "Penalty terms" = paste(penalty.control$penalty.terms, collapse = ", "),
    "Model name" = model_names,
    "Model code" = model_codes_chr,
    Fitness      = fits,
    stringsAsFactors = FALSE
  )

  best_idx <- which.min(results_table$Fitness)

  list(
    results_table = results_table,
    best_code     = candidate_codes[[best_idx]],
    best_row      = results_table[best_idx, , drop = FALSE]
  )
}


#' Screen residual error model (rv)
#'
#' Uses the current model as baseline and, if \code{rv == 3}, also evaluates
#' \code{rv = 1} and \code{rv = 2} while keeping all other fields fixed.
#' Chooses the candidate with the lowest Fitness.
#'
#' @param dat Data frame passed to \code{mod.run()}.
#' @param state List with optional \code{best_code} (named integer vector) and \code{modi}.
#' @param search.space Character scalar: \code{"ivbase"} or \code{"oralbase"}.
#' @param param_table Parameter table passed to \code{mod.run()}.
#' @param penalty.control Optional penalty control passed to \code{mod.run()}.
#' @param ... Forwarded to \code{mod.run()} (e.g., \code{filename}, controls, \code{custom_base}).
#'
#' @return List with:
#'   \itemize{
#'     \item \code{results_table}: data.frame (Step, Model name, Model code, Fitness)
#'     \item \code{best_code}: named integer vector of the best candidate
#'     \item \code{best_row}: one-row data.frame (the best candidate)
#'   }
#'
#' @examples
#' \dontrun{
#' state <- list(best_code = base_model("ivbase"), modi = 1L)
#' res_rv <- step_rv(
#'   dat = dat,
#'   state = state,
#'   search.space = "ivbase",
#'   param_table = param_table,
#'   penalty.control = penaltyControl()
#' )
#' res_rv$results_table
#' res_rv$best_code
#' }
#'
#' @export
step_rv <- function(dat,
                    state = list(),
                    search.space = "ivbase",
                    param_table = NULL,
                    penalty.control = NULL,
                    precomputed_results_file=NULL,
                    filename=NULL,
                    ...) {
  dots <- list(...)
  custom_base <-
    if (!is.null(dots$custom_base))
      dots$custom_base
  else
    NULL

  # starting code
  if (!is.null(state$best_code)) {
    current_code <- state$best_code
  } else {
    current_code <-
      base_model(search.space = search.space, custom_base = custom_base)
  }

  # baseline candidate
  candidate_codes <- list(current_code)

  candidate_codes <- lapply(1:3, function(rv_val) {
    code <- current_code
    code["rv"] <- rv_val
    code
  })

  # evaluate all candidates
  fits <- vapply(candidate_codes, function(code_vec) {
    mod.run(
      modi             = modi,
      string           = unname(code_vec),
      dat              = dat,
      search.space     = search.space,
      param_table      = param_table,
      penalty.control  = penalty.control,
      precomputed_results_file=   precomputed_results_file,
      filename=filename,
      ...
    )
  }, numeric(1))

  # names & codes
  model_names <- vapply(candidate_codes,
                        function(code)
                          CodetoMod(search.space = search.space, sel.best.code = unname(code)),
                        character(1))

  model_codes_chr <- vapply(candidate_codes, function(code) {
    paste(unname(code), collapse = ",")
  }, character(1))

  results_table <- data.frame(
    Step         = "Residual error types",
    "Penalty terms" = paste(penalty.control$penalty.terms, collapse = ", "),
    "Model name" = model_names,
    "Model code" = model_codes_chr,
    Fitness      = fits,
    stringsAsFactors = FALSE
  )

  # choose best
  best_idx <- which.min(results_table$Fitness)
  r<<-r+1
  list(
    results_table = results_table,
    best_code     = candidate_codes[[best_idx]],
    best_row      = results_table[best_idx, , drop = FALSE]
  )
}



#' Stepwise Model Building for Nonlinear Mixed-Effects Models
#'
#' The `sf.operator` function performs a stepwise model-building procedure
#' for nonlinear mixed-effects (NLME) models, including testing the number of
#' compartments, elimination type (e.g., Michaelis-Menten), inter-individual
#' variability (IIV), correlations, and residual error models.
#'
#' This function is designed to work with the `nlmixr2` and `rxode2` workflow
#' and will create intermediate model runs and results in a temporary directory.
#'
#' @param dat A dataset containing the pharmacokinetic/pharmacodynamic (PK/PD) data
#'   to be modeled.
#' @param search.space Character string specifying the structural model search space.
#'   Must be either `"ivbase"` (intravenous bolus) or `"oralbase"` (oral administration).
#' @param param_table Optional parameter table. If `NULL`, the function will generate
#'   one automatically using \code{auto_param_table()}.
#' @param penalty.control An object created by \code{penaltyControl()} specifying
#'   penalty terms to be used in model selection.
#' @param no.cores Integer specifying the number of threads for parallel execution.
#'   Defaults to \code{rxode2::getRxThreads()}.
#' @param foldername Character string for naming the output folder.
#' @param filename Character string for naming output files.
#' @param custom_base Optional custom base model to use instead of the default generated one.
#' @param dynamic_fitness Logical; if `TRUE`, penalty terms change dynamically
#'   during different steps of model building.
#' @param precomputed_results_file Optional path to a file containing precomputed
#'   stepwise search results to avoid rerunning steps.
#' @param ... Additional arguments passed to the underlying stepwise functions.
#'
#' @details
#' The procedure follows these steps:
#' \enumerate{
#'   \item **Step 1:** Select the optimal number of compartments.
#'   \item **Step 2:** Determine the elimination type (linear vs. Michaelis–Menten).
#'   \item **Step 3:** Test inclusion of inter-individual variability (IIV) terms.
#'   \item **Step 4:** Explore correlations between random effects (if any `eta.*` terms exist).
#'   \item **Step 5:** Explore different residual error models.
#' }
#'
#' If no `eta.*` parameters remain after Step 3 (sum equals zero), Step 4
#' (correlation testing) will be skipped.
#'
#' All intermediate results are stored in \code{Store.all} and the final best
#' model code is returned in a structured object of class \code{"sfOperatorResult"}.
#'
#' @return An object of class \code{"sfOperatorResult"} containing:
#' \itemize{
#'   \item \code{"Final Best Code"} – Best model code parameters.
#'   \item \code{"Final Best Model Name"} – Name of the best model.
#'   \item \code{"Stepwise Best Models"} – Summary of best models per step.
#'   \item \code{"Stepwise History"} – Detailed results from each step.
#'   \item \code{"Model Run History"} – All model runs performed.
#' }
#'
#' @seealso \code{\link{step_compartments}}, \code{\link{step_elimination}},
#'   \code{\link{step_iiv}}, \code{\link{step_correlation}}, \code{\link{step_rv}}
#'
#' @examples

#' \dontrun{
#' result<-sf.operator(dat = pheno_sd)
#' print(result)
#' }
#'
#' @export

sf.operator <- function(dat,
                        search.space = "ivbase",
                        param_table = NULL,
                        penalty.control = penaltyControl(),
                        no.cores = rxode2::getRxThreads(),
                        foldername = "test",
                        filename = "test",
                        custom_base = NULL,
                        dynamic_fitness = TRUE,
                        precomputed_results_file =NULL,
                        ...) {

  current.date<-Sys.Date()
  outputdir <-
    paste0("Step_",
           current.date,
           "-",
           foldername,
           "_",
           digest::digest(dat),
           "_temp")

  if (!dir.exists(outputdir)) {
    dir.create(outputdir, showWarnings = FALSE, recursive = TRUE)
  } else {
    message(
      sprintf(
        "Output directory '%s' already exists. Using existing directory.",
        outputdir
      )
    )
  }

  setwd(outputdir)

  # Set initial estimate
  param_table <- auto_param_table(
    dat = dat,
    param_table = param_table,
    nlmixr2autoinits = T,
    foldername = foldername
  )

  current <- base_model(search.space = search.space,
                        custom_base = custom_base)

  precomputed_cache_loaded <<- FALSE
  Store.all <<- NULL
  r <<- 1
  modi <<- 1
  #################################Step1. No. of compartment##########################
  message(crayon::blue(
    paste0(
      "Running Stepwise 1. Structural Model----------------------------------------------------"
    )
  ))

  message(crayon::blue(
    paste0(
      "Test number of compartments----------------------------------------------------"
    )
  ))

  state <- list()

   if (isTRUE(dynamic_fitness)) {
     penalty.control$penalty.terms = c("rse", "theta", "covariance")
   }

  result.steps.compartments <-   step_compartments(
    dat = dat,
    state = state,
    search.space = search.space,
    param_table = param_table,
    penalty.control = penalty.control,
    precomputed_results_file=  precomputed_results_file,
    filename=filename,
    ...
  )
  ##################### Step2. Michaelis-Menten kinetics###############
  message(crayon::blue(
    paste0(
      "Analyse elimination type----------------------------------------------------"
    )
  ))
  result.steps.MM <-   step_elimination(
    dat = dat,
    state = result.steps.compartments,
    search.space = search.space,
    param_table = param_table,
    penalty.control = penalty.control,
    filename = filename,
    ...
  )

  ####################### Step3. Random effects############################
  message(crayon::blue(
    paste0(
      "Test IIV on parameters----------------------------------------------------"
    )
  ))

  if (isTRUE(dynamic_fitness)) {
    penalty.control$penalty.terms = c("rse", "theta", "covariance", "shrinkage", "omega")
  }

  result.steps.iiv <-   step_iiv(
    dat = dat,
    state = result.steps.MM,
    search.space = search.space,
    param_table = param_table,
    penalty.control = penalty.control,
    filename = filename,
    ...
  )

  ######################## Step 3. Explore correlation. ##########################
  message(crayon::blue(
    paste0(
      "Test Correlation between parameters----------------------------------------------------"
    )
  ))

 if (isTRUE(dynamic_fitness)) {
   penalty.control$penalty.terms = c(
     "rse",
     "theta",
     "covariance",
     "shrinkage",
     "omega",
     "correlation"
   )
  }

  eta_core_sum <-
    sum(result.steps.iiv$best_code[c("eta.vc", "eta.vp", "eta.vp2", "eta.q", "eta.q2")], na.rm = TRUE)
   mm <- result.steps.iiv$best_code["mm"]

  result.steps.corr <- NULL
  if ((mm == 0 && eta_core_sum > 0) ||
      (mm == 1 && (result.steps.iiv$best_code["eta.km"] == 1 || eta_core_sum > 1))) {

    result.steps.corr <- step_correlation(
      dat = dat,
      state = result.steps.iiv,
      search.space = search.space,
      param_table = param_table,
      penalty.control = penalty.control,
      filename = filename,
      ...
    )
  }
  ######################## Step 4. Explore RV model ##############################

  message(crayon::blue(
    paste0(
      "Explore types of residual errors----------------------------------------------------"
    )
  ))

  if (!is.null( result.steps.corr)){
  result.steps.rv <-   step_rv(
    dat = dat,
    state = result.steps.corr,
    search.space = search.space,
    param_table = param_table,
    penalty.control = penalty.control,
    ...
  )
 }

  if (is.null( result.steps.corr)){
    result.steps.rv <-   step_rv(
      dat = dat,
      state = result.steps.iiv,
      search.space = search.space,
      param_table = param_table,
      penalty.control = penalty.control,
      ...
    )
  }

  out <- new.env(parent = emptyenv())
  class(out) <- "sfOperatorResult"
  latest_round <- subset(Store.all, round.num == max(Store.all$round.num, na.rm = TRUE))
  best_row <- latest_round[which.min(latest_round$fitness), ]

  if (search.space == "ivbase") {
    cols_to_extract <- c(
      "no.cmpt", "eta.vmax", "eta.km", "eta.cl", "eta.vc",
      "eta.vp", "eta.vp2", "eta.q", "eta.q2",
      "mm", "mcorr", "rv"
    )
  } else {
    cols_to_extract <- c(
      "no.cmpt", "eta.vmax", "eta.km", "eta.cl", "eta.vc",
      "eta.vp", "eta.vp2", "eta.q", "eta.q2", "eta.ka",
      "mm", "mcorr", "rv"
    )
  }

  # Assign the selected columns from the best row to the output
  out[["Final Best Code"]] <- best_row[, cols_to_extract, drop = FALSE]

  out[["Final Best Model Name"]] <-
    result.steps.rv$best_row$Model.name
  out[["Stepwise Best Models"]] <- data.frame(
    rbind(
      result.steps.compartments$best_row,
      result.steps.MM$best_row,
      result.steps.iiv$best_row,
      result.steps.corr$best_row,
      result.steps.rv$best_row
    ),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  out[["Stepwise History"]]     <- list(
    "Step 1 Compartments" = result.steps.compartments,
    "Step 2 Elimination"  = result.steps.MM,
    "Step 3 IIV"          = result.steps.iiv,
    "Step 4 Correlation"  = result.steps.corr,
    "Step 5 RV"           = result.steps.rv
  )

  out[["Model Run History"]] <-
    as.data.frame(Store.all, stringsAsFactors = FALSE)

  return(out)
}


#' Print Method for sfOperatorResult Objects
#'
#' Defines a custom print method for objects of class
#' \code{sfOperatorResult}. It displays the best model code, best model name,
#' and the stepwise selection history with formatted colors using the
#' \pkg{crayon} package.
#'
#' @param x An object of class \code{sfOperatorResult}, typically returned by
#'   a model selection procedure.
#' @param ... Further arguments passed to or from other methods (currently unused).
#'
#' @return Invisibly returns \code{x}, after printing its contents to the console.
#'
#' @details
#' This method prints:
#' \itemize{
#'   \item \strong{Best Model Code}: The final chosen model's code.
#'   \item \strong{Best Model Name}: The final chosen model's name.
#'   \item \strong{Stepwise Selection History}: A record of the best models
#'         found at each step of a stepwise selection process.
#' }
#'
#' @examples
#' \dontrun{
#' # Assuming `res` is an object of class "sfOperatorResult"
#' print(res)
#' }
#'
#' @seealso \code{\link[crayon]{crayon}} for text formatting utilities.
#'
#' @export
print.sfOperatorResult <- function(x, ...) {
  cat(crayon::green$bold("\n=== Best Model Code ===\n"))
  print(x$`Final Best Code`)
  cat(crayon::green$bold("\n=== Best Model Name ===\n"))
  cat(x$`Final Best Model Name`, "\n")
  cat("\n=== Stepwise Selection History ===\n")
  print(x$`Stepwise Best Models`)
}


