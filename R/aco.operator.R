
#' Initialize Node List for ACO Search Space
#'
#' Creates the initial `node.list` data frame for the Ant Colony Optimization (ACO) search,
#' defining nodes, their groups, initial pheromone values, and selection probabilities.
#' Supports two search spaces: intravenous ("ivbase") and oral ("oralbase").
#'
#' @param search.space Character. Type of search space: `"ivbase"` (intravenous) or `"oralbase"` (oral, includes `eta.ka` nodes).
#' @param initial.phi Numeric. Initial pheromone value assigned to each node.
#'
#' @return A data frame representing the initialized node list, including:
#'   \itemize{
#'     \item `travel` — Travel counter (initialized to 0).
#'     \item `node.no` — Node index number.
#'     \item `local.node.no` — Index of the option within the node group.
#'     \item `node.names` — Node name (e.g., `"eta.vp.no"`, `"eta.vp.yes"`).
#'     \item `node.group` — Group ID for nodes that share a pheromone distribution.
#'     \item `phi` — Initial pheromone value.
#'     \item `delta_phi` — Change in pheromone (initialized to 0).
#'     \item `p` — Initial selection probability for each node option.
#'   }
#' @examples
#' node.list.iv <- initNodeList(search.space = "ivbase", initial.phi = 1)
#' node.list.oral <- initNodeList(search.space = "oralbase", initial.phi = 1)
#'
#' @export
initNodeList <- function(search.space,
                         initial.phi) {
  if (search.space == "ivbase") {
    no.nodes <- 22
    node.list.all <- data.frame(
      travel = 0,
      node.no = seq(1, no.nodes, 1),
      local.node.no = c(
        seq(1, 3, 1),
        # compartments: 1Cmpt / 2Cmpt / 3Cmpt
        seq(0, 1, 1),
        # vp2
        seq(0, 1, 1),
        # q2
        seq(0, 1, 1),
        # vp
        seq(0, 1, 1),
        # q
        seq(0, 1, 1),
        # vc
        seq(0, 1, 1),
        # mm
        seq(0, 1, 1),
        # km
        seq(0, 1, 1),
        # mcorr
        seq(1, 3, 1)   # residual error: add / prop / comb
      ),
      node.names = c(
        "1Cmpt",
        "2Cmpt",
        "3Cmpt",
        "eta.vp2.no",
        "eta.vp2.yes",
        "eta.q2.no",
        "eta.q2.yes",
        "eta.vp.no",
        "eta.vp.yes",
        "eta.q.no",
        "eta.q.yes",
        "eta.vc.no",
        "eta.vc.yes",
        "mm.no",
        "mm.yes",
        "eta.km.no",
        "eta.km.yes",
        "mcorr.no",
        "mcorr.yes",
        "add",
        "prop",
        "comb"
      ),
      node.group = c(rep(1, 3),
                     sort(rep(seq(
                       2, 9, 1
                     ), 2)),
                     rep(10, 3)),
      phi = rep(initial.phi, no.nodes),
      delta_phi = rep(0, no.nodes),
      p = c(rep(round(1 / 3, 3), 3),
            rep(0.5, 16),
            rep(round(1 / 3, 3), 3))
    )

  } else if (search.space == "oralbase") {
    no.nodes <- 24
    node.list.all <- data.frame(
      travel = 0,
      node.no = seq(1, no.nodes, 1),
      local.node.no = c(
        seq(1, 3, 1),
        # compartments: 1Cmpt / 2Cmpt / 3Cmpt
        seq(0, 1, 1),
        # vp2
        seq(0, 1, 1),
        # q2
        seq(0, 1, 1),
        # vp
        seq(0, 1, 1),
        # q
        seq(0, 1, 1),
        # vc
        seq(0, 1, 1),
        # ka (oral-specific absorption rate constant)
        seq(0, 1, 1),
        # mm
        seq(0, 1, 1),
        # km
        seq(0, 1, 1),
        # mcorr
        seq(1, 3, 1)   # residual error: add / prop / comb
      ),
      node.names = c(
        "1Cmpt",
        "2Cmpt",
        "3Cmpt",
        "eta.vp2.no",
        "eta.vp2.yes",
        "eta.q2.no",
        "eta.q2.yes",
        "eta.vp.no",
        "eta.vp.yes",
        "eta.q.no",
        "eta.q.yes",
        "eta.vc.no",
        "eta.vc.yes",
        "eta.ka.no",
        "eta.ka.yes",
        # oral-specific
        "mm.no",
        "mm.yes",
        "eta.km.no",
        "eta.km.yes",
        "mcorr.no",
        "mcorr.yes",
        "add",
        "prop",
        "comb"
      ),
      node.group = c(rep(1, 3),
                     sort(rep(seq(
                       2, 10, 1
                     ), 2)),
                     rep(11, 3)),
      phi = rep(initial.phi, no.nodes),
      delta_phi = rep(0, no.nodes),
      p = c(rep(round(1 / 3, 3), 3),
            rep(0.5, (no.nodes - 6)),
            rep(round(1 / 3, 3), 3))
    )

  } else {
    stop("Unknown search.space type: must be 'ivbase' or 'oralbase'")
  }

  return(node.list.all)
}



#' Control Parameters for Ant Colony Optimization
#'
#' Creates a list of control settings for the \code{\link{aco.operator}} function.
#' These settings define the number of ants, iterations, pheromone update rules,
#' probability bounding, and elitism.
#'
#' @param no.ants Integer. Number of ants per iteration.
#' @param max.iter Integer. Maximum number of iterations.
#' @param rho Numeric in [0,1]. Pheromone evaporation rate per iteration.
#' @param initial.phi Numeric. Initial pheromone level for all nodes.
#' @param lower.limit.phi Numeric. Lower bound for pheromone level.
#' @param upper.limit.phi Numeric. Upper bound for pheromone level.
#' @param alpha.value Numeric. Influence factor of pheromone level on path selection.
#' @param elitism.percentage Numeric in [0,1]. Fraction of top ants preserved each iteration.
#' @param prob.floor Numeric in (0, 0.5). Minimum selection probability assigned to
#' each node option (probability floor constraint). After applying the floor, the
#' remaining probability mass is proportionally rescaled among the other options
#' to preserve relative likelihoods.
#' @param sig.diff Numeric. Minimum significant fitness difference for ranking.
#'
#' @details
#' The \code{prob.floor} parameter implements a *probability floor constraint with proportional
#' rescaling* (also known as bounded probability adjustment). This prevents any decision
#' option from having a probability lower than \code{prob.floor}, improving exploration
#' stability in cases where pheromone values become highly imbalanced.
#'
#' @return List of ACO hyperparameters.
#' @export

acoControl <- function(no.ants = 10,
                       max.iter = 10,
                       rho = 0.2,
                       initial.phi = 1,
                       lower.limit.phi = 1,
                       upper.limit.phi = Inf,
                       alpha.value = 1,
                       elitism.percentage = 1 / 3,
                       prob.floor = 0.2,
                       sig.diff = 1) {
  list(
    no.ants = no.ants,
    max.iter = max.iter,
    rho = rho,
    initial.phi = initial.phi,
    lower.limit.phi = lower.limit.phi,
    upper.limit.phi = upper.limit.phi,
    alpha.value = alpha.value,
    elitism.percentage = elitism.percentage,
    prob.floor = prob.floor,
    sig.diff = sig.diff
  )
}


#' Ant Colony Optimization (ACO) Operator for Model Selection
#'
#' Implements an Ant Colony Optimization (ACO) algorithm to explore
#' model space (e.g., pharmacometric structural models) and identify the best-performing
#' model given observed data and parameter constraints.
#'
#' The ACO approach uses a colony of "ants" to stochastically sample models,
#' evaluate their fitness, and update pheromone trails that guide future searches.
#' This iterative process balances exploration of new models with exploitation
#' of promising candidates.
#'
#' @param dat A dataset (typically pharmacokinetic/pharmacodynamic data) used
#'   for model estimation and evaluation.
#' @param param_table Optional parameter table specifying initial estimates,
#'   bounds, and parameter definitions. If `NULL`, it will be automatically
#'   generated via \code{auto_param_table()}.
#' @param search.space Character string specifying the search space type.
#'   Options are:
#'   \itemize{
#'     \item `"ivbase"` — base IV model search space.
#'     \item `"oralbase"` — base oral model search space.
#'   }
#' @param no.cores Integer. Number of CPU cores to use for parallelization.
#'   Defaults to \code{rxode2::getRxThreads()}.
#' @param foldername Character string. Base folder name for output storage.
#' @param filename Character string. Base file name for result storage.
#' @param aco.control A list of ACO control parameters, typically generated by
#'   \code{acoControl()}. Includes:
#'   \itemize{
#'     \item no.ants — number of ants per iteration.
#'     \item max.iter — maximum number of iterations.
#'     \item rho — pheromone evaporation rate.
#'     \item initial.phi — initial pheromone level.
#'     \item lower.limit.phi, upper.limit.phi — bounds for pheromone levels.
#'     \item alpha.value — pheromone weight exponent.
#'     \item elitism.percentage — proportion of best solutions preserved each iteration.
#'     \item prob.floor — minimum sampling probability.
#'     \item sig.diff — threshold for significant fitness difference.
#'   }
#' @param penalty.control A list of penalty function settings, typically from
#'   \code{penaltyControl()}, to handle invalid or poor-fitting models.
#' @param precomputed_results_file Optional. A file path to precomputed model results
#'   to avoid redundant evaluations.
#' @param seed.no Integer. Random seed for reproducibility.
#' @param ... Additional arguments passed to the underlying model fitting function
#'   (e.g., \code{mod.run()}).
#'
#' @return An object of class \code{"acoOperatorResult"}, containing:
#' \itemize{
#'   \item \code{$`Final Selected Code`} — Vector representation of the best model.
#'   \item \code{$`Final Selected Model Name`} — Human-readable name of the selected model.
#'   \item \code{$`Model Run History`} — Data frame of model runs across iterations.
#'   \item \code{$`Node Run History`} — History of pheromone probabilities for each iteration.
#' }
#'
#' @examples
#' \dontrun{
#' # Example usage with phenotype dataset
#' outs <- aco.operator(
#'   dat = pheno_sd,
#'   param_table = NULL,
#'   search.space = "ivbase",
#'   saem.control = saemControl(
#'     seed = 1234,
#'     nBurn = 15,
#'     nEm   = 15,
#'     logLik = TRUE
#'   )
#' )
#'
#' # Print summary
#' print(outs)
#' }
#'
#' @seealso \code{\link{acoControl}}, \code{\link{penaltyControl}}, \code{\link{saemControl}}
#'
#' @export

aco.operator <- function(dat,
                         param_table = NULL,
                         search.space = "ivbase",
                         no.cores = rxode2::getRxThreads(),
                         foldername = "test",
                         filename = "test",
                         aco.control = acoControl(),
                         penalty.control = penaltyControl(),
                         precomputed_results_file = NULL,
                         seed.no = 1234,
                         ...) {
  # --- Extract ACO control parameters ---
  no.ants <- aco.control$no.ants
  max.iter <- aco.control$max.iter
  rho <- aco.control$rho
  initial.phi <- aco.control$initial.phi
  lower.limit.phi <- aco.control$lower.limit.phi
  upper.limit.phi <- aco.control$upper.limit.phi
  alpha.value <- aco.control$alpha.value
  elitism.percentage <- aco.control$elitism.percentage
  prob.floor <- aco.control$prob.floor
  sig.diff <- aco.control$sig.diff
  # --- Setup ---
  set.seed(seed.no)
  current.date <- Sys.Date()

  # Create temporary output directory
  outputdir <- paste0("ACO_",
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

  # Search Space Definiation
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

  # --- Iterative Tabu Search --
  progressr::handlers(
    progressr::handler_progress(
      format = paste0(
        crayon::cyan("ACO Search "),
        crayon::yellow("[:bar]"),
        crayon::green(" :percent "),
        crayon::blue(" (iteration :current/:total)")
      ),
      width = 80
    )
  )

  with_progress({
    p <- progressr::progressor(steps = max.iter)
    ##########################Start###############################
    #  Subsequent iterations ---
    for (aco.iter in 1:max.iter) {
      if (aco.iter == 1) {
        node.list.0 <- initNodeList(search.space = search.space,
                                    initial.phi = initial.phi)
        node.list.history <- node.list.0
        cycle.all.list <- list()

        # Create initial ant population
        initial.ants <- createAnts(
          search.space = "ivbase",
          no.ants = no.ants,
          initialize = TRUE,
          node.list = node.list.0
        )

        initial.ants <- t(vapply(seq_len(ncol(initial.ants)),  #
                                 function(i) {
                                   validateModels(
                                     string = pmax(unname(initial.ants[, i]), 0),
                                     search.space = "ivbase",
                                     code.source = "ACO"
                                   )
                                 },
                                 numeric(nrow(initial.ants))))

        colnames(initial.ants) <- bit.names
        cycle.all.list[[1]] <- initial.ants

        # Model running and fitness evaluation
        data.ants <- as.data.frame(initial.ants)

        data.ants$fitness <- vapply(seq_len(nrow(data.ants)),
                                    function(k) {
                                      string_vec <- as.vector(initial.ants[k, ])
                                      result <- try(mod.run(
                                        r                = aco.iter,
                                        dat              = dat,
                                        search.space     = search.space,
                                        string           = string_vec,
                                        param_table      = param_table,
                                        penalty.control  = penalty.control,
                                        precomputed_results_file = precomputed_results_file,
                                        filename         = filename,
                                        ...
                                      ),
                                      silent = TRUE)
                                      if (is.numeric(result) &&
                                          length(result) == 1)
                                        result
                                      else
                                        NA_real_
                                    },
                                    numeric(1))

        data.ants$round <- aco.iter
        data.ants$mod.no <- seq_len(nrow(data.ants))
        data.ants <-
          data.ants[, c("round", "mod.no", setdiff(names(data.ants), c("round", "mod.no")))]

        fitness_history <- data.frame()
        data.ants$round <-
          aco.iter  # Track which generation each record belongs to
        fitness_history <- rbind(fitness_history, data.ants)

        #  Update pheromone levels and probabilities ---
        node.list.1 <- phi.calculate(
          r = aco.iter,
          search.space = search.space,
          fitness_history = fitness_history,
          node.list.history = node.list.history,
          alpha.value = alpha.value,
          rho = rho,
          sig.diff = sig.diff,
          lower.limit.phi = lower.limit.phi,
          upper.limit.phi = upper.limit.phi
        )

        node.list.s <- p.calculation(node.list = node.list.1,
                                     prob.floor = prob.floor)

        node.list.history <- rbind(node.list.history, node.list.s)
      } else{
        # Identify current best model (elitism)
        bestmodel <-
          fitness_history[fitness_history$fitness == min(fitness_history$fitness), ]
        bestmodelcode <- bestmodel[1, bit.names]

        # Generate next generation of ants
        cycle.all <- createAnts(
          search.space = search.space,
          no.ants = no.ants,
          initialize = FALSE,
          node.list = node.list.s
        )

        # Apply elitism: preserve top-performing ants
        no.elitism <- max(round(no.ants * elitism.percentage, 0), 1)
        for (loop.no.elitism in 1:no.elitism) {
          cycle.all[, (no.ants - loop.no.elitism + 1)] <-
            as.numeric(bestmodelcode)
        }
        cycle.all.list[[aco.iter]] <- cycle.all

        # Evaluate all ants in current iteration
        data.ants <- as.data.frame(t(vapply(seq_len(ncol(cycle.all)),
                                            function(i) {
                                              validateModels(
                                                string = pmax(unname(cycle.all[, i]), 0),
                                                search.space = search.space,
                                                code.source = "ACO"
                                              )
                                            },
                                            numeric(nrow(cycle.all)))))
        colnames(data.ants) <- bit.names

        data.ants$fitness <- vapply(seq_len(nrow(data.ants)),
                                    function(k) {
                                      string_vec <- as.vector(as.numeric(data.ants[k, ]))
                                      result <- try(mod.run(
                                        r = aco.iter,
                                        dat = dat,
                                        search.space = search.space,
                                        string = string_vec,
                                        param_table = param_table,
                                        penalty.control = penalty.control,
                                        precomputed_results_file = precomputed_results_file,
                                        filename = filename,
                                        ...
                                      ),
                                      silent = TRUE)
                                      if (is.numeric(result) &&
                                          length(result) == 1)
                                        result
                                      else
                                        NA_real_
                                    },
                                    numeric(1))

        # Add round and model IDs
        data.ants$round <- aco.iter
        data.ants$mod.no <- seq_len(nrow(data.ants))
        data.ants <-
          data.ants[, c("round", "mod.no", setdiff(names(data.ants), c("round", "mod.no")))]

        # Append to fitness history
        fitness_history <-
          rbind(fitness_history[, setdiff(names(fitness_history), "allrank")],
                data.ants)

        # Update pheromone trails
        node.list.s <- phi.calculate(
          r = aco.iter,
          search.space = search.space,
          fitness_history = fitness_history,
          node.list.history = node.list.history,
          alpha.value = alpha.value,
          rho = rho,
          sig.diff = sig.diff,
          lower.limit.phi = lower.limit.phi,
          upper.limit.phi = upper.limit.phi
        )

        # Update probabilities
        node.list.s <- p.calculation(node.list = node.list.s,
                                     prob.floor = prob.floor)

        # Extend node history
        node.list.history <- rbind(node.list.history, node.list.s)
      }
    }
  })
  # ----------------------------
  # Final output (ACO)
  # ----------------------------
  names(bestmodelcode) <- bit.names
  best_model_name <- CodetoMod(sel.best.code = bestmodelcode,
                               search.space  = search.space)

  out <- new.env(parent = emptyenv())
  class(out) <- "acoOperatorResult"
  out[["Final Selected Code"]] <- bestmodelcode
  out[["Final Selected Model Name"]] <- best_model_name
  out[["Model Run History"]] <-
    as.data.frame(Store.all, stringsAsFactors = FALSE)
  out[["Node Run History"]] <- node.list.history

  on.exit({
    rm(modi, r, Store.all, precomputed_cache_loaded, envir = .GlobalEnv)
  }, add = TRUE)

  return(out)
}

#' Print Method for ACO Operator Results
#'
#' This function defines the \code{print} method for objects of class
#' \code{"acoOperatorResult"} produced by \code{\link{aco.operator}}.
#' It displays the final selected model code and the corresponding model name
#' in a formatted, colorized output.
#'
#' @param x An object of class \code{"acoOperatorResult"}, typically the output
#'   of \code{\link{aco.operator}}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns the input object \code{x}.
#'
#' @examples
#' \dontrun{
#' # Run ACO model selection
#' outs <- aco.operator(dat = pheno_sd)
#'
#' # Print results in formatted style
#' print(outs)
#' }
#'
#' @seealso \code{\link{aco.operator}}
#'
#' @export

print.acoOperatorResult <- function(x, ...) {
  # Print final selected model code
  cat(crayon::blue$bold("\n=== Final Selected Model Code (ACO) ===\n"))
  print(x$`Final Selected Code`)

  # Print final selected model name
  cat(crayon::blue$bold("\n=== Final Selected Model Name (ACO) ===\n"))
  cat(x$`Final Selected Model Name`, "\n")

  invisible(x)
}
