#' Initialize node list for ACO search space
#'
#' Construct the initial edge list used in model structure search
#' based on ant colony optimization.
#'
#' @param search.space Character, one of "ivbase" or "oralbase".
#'   Default is "ivbase".
#' @param phi0 A non-negative numeric value. Initial pheromone value assigned
#' to all nodes at the start of the search. Defaults to 2.
#' @return A data.frame in which each row represents an edge in the ACO
#'   path-construction graph, with the following columns:
#'   \describe{
#'     \item{travel}{Integer. Travel counter associated with the edge,
#'       initialized to zero.}
#'     \item{node.no}{Integer. Decision node identifier corresponding to
#'       a model feature.}
#'     \item{node.name}{Character. Semantic label of the decision node.}
#'     \item{edge.no}{Integer. Global edge index.}
#'     \item{local.edge.no}{Integer. Index of the edge within the
#'       corresponding decision node.}
#'     \item{edge.name}{Character. Semantic label of the edge (model
#'       component choice).}
#'     \item{phi}{Numeric. Initial pheromone value associated with the
#'       edge.}
#'     \item{delta_phi}{Numeric. Change in pheromone level, initialized
#'       to zero.}
#'     \item{p}{Numeric. Initial selection probability of the edge.}
#'   }
#'
#' @author Zhonghui Huang
#'
#' @examples
#' initNodeList(search.space = "ivbase", phi0 = 1)
#' initNodeList(search.space = "oralbase", phi0 = 1)
#'
#' @export

initNodeList <- function(search.space, phi0) {
  if (!search.space %in% c("ivbase", "oralbase")) {
    stop("search.space must be 'ivbase' or 'oralbase'")
  }
  # candidate nodes
  candidate_nodes <- list(
    compartments = list(
      applicable = c("ivbase", "oralbase"),
      nodes = data.frame(
        edge.name    = c("1Cmpt", "2Cmpt", "3Cmpt"),
        local.edge.no = 1:3,
        node.name      = "compartments",
        stringsAsFactors = FALSE
      )
    ),

    eta_vp2 = list(
      applicable = c("ivbase", "oralbase"),
      nodes = data.frame(
        edge.name    = c("eta.vp2.no", "eta.vp2.yes"),
        local.edge.no = 0:1,
        node.name      = "eta_vp2",
        stringsAsFactors = FALSE
      )
    ),

    eta_q2 = list(
      applicable = c("ivbase", "oralbase"),
      nodes = data.frame(
        edge.name    = c("eta.q2.no", "eta.q2.yes"),
        local.edge.no = 0:1,
        node.name      = "eta_q2",
        stringsAsFactors = FALSE
      )
    ),

    eta_vp = list(
      applicable = c("ivbase", "oralbase"),
      nodes = data.frame(
        edge.name    = c("eta.vp.no", "eta.vp.yes"),
        local.edge.no = 0:1,
        node.name      = "eta_vp",
        stringsAsFactors = FALSE
      )
    ),

    eta_q = list(
      applicable = c("ivbase", "oralbase"),
      nodes = data.frame(
        edge.name    = c("eta.q.no", "eta.q.yes"),
        local.edge.no = 0:1,
        node.name      = "eta_q",
        stringsAsFactors = FALSE
      )
    ),

    eta_vc = list(
      applicable = c("ivbase", "oralbase"),
      nodes = data.frame(
        edge.name    = c("eta.vc.no", "eta.vc.yes"),
        local.edge.no = 0:1,
        node.name      = "eta_vc",
        stringsAsFactors = FALSE
      )
    ),

    eta_ka = list(
      applicable = "oralbase",
      nodes = data.frame(
        edge.name    = c("eta.ka.no", "eta.ka.yes"),
        local.edge.no = 0:1,
        node.name      = "eta_ka",
        stringsAsFactors = FALSE
      )
    ),

    mm = list(
      applicable = c("ivbase", "oralbase"),
      nodes = data.frame(
        edge.name    = c("mm.no", "mm.yes"),
        local.edge.no = 0:1,
        node.name      = "mm",
        stringsAsFactors = FALSE
      )
    ),

    eta_km = list(
      applicable = c("ivbase", "oralbase"),
      nodes = data.frame(
        edge.name    = c("eta.km.no", "eta.km.yes"),
        local.edge.no = 0:1,
        node.name      = "eta_km",
        stringsAsFactors = FALSE
      )
    ),

    mcorr = list(
      applicable = c("ivbase", "oralbase"),
      nodes = data.frame(
        edge.name    = c("mcorr.no", "mcorr.yes"),
        local.edge.no = 0:1,
        node.name      = "mcorr",
        stringsAsFactors = FALSE
      )
    ),

    residual = list(
      applicable = c("ivbase", "oralbase"),
      nodes = data.frame(
        edge.name    = c("add", "prop", "comb"),
        local.edge.no = 1:3,
        node.name      = "residual",
        stringsAsFactors = FALSE
      )
    )
  )

  offset <- if (search.space == "ivbase") 0 else 1

  group_map <- c(
    compartments = 1,
    eta_vp2      = 2,
    eta_q2       = 3,
    eta_vp       = 4,
    eta_q        = 5,
    eta_vc       = 6,
    if (search.space == "oralbase") c(eta_ka = 7),
    mm           = 7 + offset,
    eta_km       = 8 + offset,
    mcorr        = 9 + offset,
    residual     = 10 + offset
  )

  node_tables <- lapply(candidate_nodes, function(x) {
    if (search.space %in% x$applicable) {
      df <- x$nodes
      df$node.no <- unname(group_map[df$node.name[1]])
      df
    } else {
      NULL
    }
  })

  node_tables <- Filter(Negate(is.null), node_tables)
  node.list <- do.call(rbind, node_tables)

  node.list <- node.list[order(node.list$node.no,
                               node.list$local.edge.no), ]

  node.list$travel    <- 0
  node.list$edge.no   <- seq_len(nrow(node.list))
  node.list$phi       <- phi0
  node.list$delta_phi <- 0

  node.list$p <- NA_real_
  node.list$p[node.list$node.name == "compartments"] <-
    round(1 / 3, 3)
  node.list$p[node.list$node.name == "residual"]     <-
    round(1 / 3, 3)
  node.list$p[is.na(node.list$p)]                   <- 0.5

  node.list <- node.list[, c(
    "travel",
    "node.no",
    "node.name",
    "edge.no",
    "local.edge.no",
    "edge.name",
    "phi",
    "delta_phi",
    "p"
  )]

  rownames(node.list) <- NULL
  return(node.list)
}


#' Create ant population for ACO
#'
#' Generate a population of ants (candidate models) for use in an
#' ant colony optimization algorithm for pharmacometric model search.
#'
#' @param search.space Character, one of "ivbase" or "oralbase".
#'   Default is "ivbase".
#' @param nants Integer. Number of ants (candidate solutions) generated
#'   at each iteration. Defaults to 15.
#' @param init Logical. If TRUE, a subset of ants is initialized as fixed base
#'   models and the remaining ants are generated by probabilistic sampling.
#'   If FALSE, all ants are generated by probabilistic sampling.
#' @param nodeslst Data frame containing pheromone information used for
#'   probabilistic sampling. It must include node identifiers and associated
#'   sampling probabilities. This argument is required whenever random ants
#'   are generated.
#' @param custom_config Optional named list defining a custom parameter structure.
#'   If provided, the parameter names are taken from the names of this list.
#'   If NULL, a default parameter structure is used based on the selected
#'   search space.
#' @param fixed Optional list specifying fixed ants for initialization.
#'   The list may contain the following elements:
#'   \itemize{
#'     \item n: number of fixed ants.
#'     \item mat: optional matrix specifying fixed ant encodings, with
#'       parameters in rows and ants in columns.
#'   }
#'
#' @details
#' Each ant is represented as a column vector encoding discrete structural
#' model decisions, including the number of compartments, inclusion of random
#' effects, Michaelis--Menten elimination, correlation structures, and residual
#' error models. The set of parameters included in the encoding depends on the
#' selected search space.
#'
#' Ants are generated using a combination of fixed initialization and
#' pheromone-guided probabilistic sampling. When fixed initialization is
#' enabled, a subset of ants corresponds to predefined base models, such as
#' one- to three-compartment structures with different residual error models.
#' The remaining ants are sampled according to probability distributions
#' derived from pheromone weights stored in the node list.
#'
#' Structural dependencies between parameters are enforced during generation.
#' For example, parameters associated with peripheral compartments are only
#' active when the number of compartments is sufficient, and parameters related
#' to Michaelis--Menten elimination are only sampled when the corresponding
#' mechanism is selected. Parameters that are not applicable for a given
#' structure are encoded with a value of -1.
#'
#' @return
#' A numeric matrix in which rows correspond to model parameters and columns
#' correspond to individual ants. Column names identify ants sequentially.
#'
#' Parameter values are encoded as integers. Binary indicators represent
#' exclusion or inclusion of model components, categorical values represent
#' multi-level structural choices, and the value -1 indicates that a parameter
#' is not applicable for the given model structure.
#'
#' @author Zhonghui Huang
#'
#' @examples
#' # Example 1: Use defaults (6 fixed base models)
#' nodes <- initNodeList(search.space = "ivbase", phi0 = 1)
#' createAnts(
#'   search.space = "ivbase",
#'   nants = 15,
#'   init = TRUE,
#'   nodeslst = nodes
#' )
#'
#' # Example 2: Custom number of fixed ants
#' createAnts(
#'   search.space = "ivbase",
#'   nants = 20,
#'   init = TRUE,
#'   nodeslst = nodes,
#'   fixed = list(n = 10)  # 10 fixed, mat = NULL (auto-generate)
#' )
#'
#' # Example 3: Custom fixed models
#' my_models <- matrix(
#'   c(1, -1, 0, -1, -1, -1, -1, 0, 0, 1,
#'     2, -1, 1,  0, -1,  0, -1, 0, 0, 2,
#'     3, -1, 1,  1,  0,  1,  0, 1, 1, 3),
#'   nrow = 10, ncol = 3
#' )
#' rownames(my_models) <- c("no.cmpt", "eta.km", "eta.vc", "eta.vp",
#'                          "eta.vp2", "eta.q", "eta.q2", "mm", "mcorr", "rv")
#'
#' createAnts(
#'   search.space = "ivbase",
#'   nants = 10,
#'   init = TRUE,
#'   nodeslst = nodes,
#'   fixed = list(n = 3, mat = my_models)
#' )
#' @seealso \code{\link{initNodeList}}, \code{\link{aco.operator}}
#'
#' @export

createAnts <- function(search.space = "ivbase",
                       nants = 15,
                       init = FALSE,
                       nodeslst = NULL,
                       custom_config = NULL,
                       fixed = NULL) {

  if (is.null(custom_config)) {
    space_cfg <- spaceConfig(search.space)
    params <- space_cfg$params
  } else {
    params <- names(custom_config)
  }

  np <- length(params)

  if (init) {
    # Setup fixed configuration
    fixed_config <- list(n = 6,
                         mat = NULL)

    if (!is.null(fixed)) {
      if ("n" %in% names(fixed)) {
        fixed_config$n <- fixed$n
      }
      if ("mat" %in% names(fixed)) {
        fixed_config$mat <- fixed$mat
      }
    }
    nf <- fixed_config$n
    if (nants < nf) {
      warning("nants (",
              nants,
              ") < n_fixed (",
              nf,
              "). Only ",
              nants,
              " fixed ants.")
      nf <- nants
      nr <- 0
    } else {
      nr <- nants - nf
    }

    # Create fixed ants matrix
    if (!is.null(fixed_config$mat)) {
      # Custom matrix
      fx <- fixed_config$mat
      if (nrow(fx) != np) {
        stop("fixed$mat must have ", np, " rows")
      }
      if (ncol(fx) != nf) {
        stop("fixed$mat must have ", nf, " columns")
      }
      rownames(fx) <- params
      colnames(fx) <- paste0("ant", 1:nf)

    } else {
      # Auto-generate default fixed models
      fx <- matrix(NA, nrow = np, ncol = nf)
      rownames(fx) <- params
      colnames(fx) <- paste0("ant", 1:nf)

      # Pattern: cycle through 1/2/3 compartments with add/prop errors
      cmpt <- rep(1:3, each = 2, length.out = nf)
      rv <- rep(c(1, 2), length.out = nf)

      fx["no.cmpt", ] <- cmpt
      fx["rv", ] <- rv
      fx["eta.vc", ] <- 0
      fx["mm", ] <- 0
      fx["mcorr", ] <- 0
      fx["eta.km", ] <- -1

      for (i in 1:nf) {
        fx["eta.vp", i] <- if (cmpt[i] >= 2)
          0
        else
          - 1
        fx["eta.q", i] <- if (cmpt[i] >= 2)
          0
        else
          - 1
        fx["eta.vp2", i] <- if (cmpt[i] == 3)
          0
        else
          - 1
        fx["eta.q2", i] <- if (cmpt[i] == 3)
          0
        else
          - 1
      }

      if ("eta.ka" %in% params) {
        fx["eta.ka", ] <- -1
      }
    }

  } else {
    # No fixed ants
    nf <- 0
    nr <- nants
    fx <- NULL
  }

  if (nr > 0) {
    if (is.null(nodeslst)) {
      stop("nodeslst required for random ant generation")
    }
    # Create random ants matrix
    rx <- matrix(NA, nrow = np, ncol = nr)
    rownames(rx) <- params

    idx_cmpt <- grep("Cmpt$", nodeslst$edge.name)
    rx["no.cmpt", ] <-
      sample(1:3, nr, prob = nodeslst$p[idx_cmpt], replace = TRUE)

    idx_vc <- grep("^eta\\.vc\\.", nodeslst$edge.name)
    rx["eta.vc", ] <-
      sample(0:1, nr, prob = nodeslst$p[idx_vc], replace = TRUE)

    idx_mm <- grep("^mm\\.", nodeslst$edge.name)
    rx["mm", ] <-
      sample(0:1, nr, prob = nodeslst$p[idx_mm], replace = TRUE)

    idx_mc <- grep("^mcorr\\.", nodeslst$edge.name)
    rx["mcorr", ] <-
      sample(0:1, nr, prob = nodeslst$p[idx_mc], replace = TRUE)

    idx_rv <- grep("^(add|prop|comb)$", nodeslst$edge.name)
    rx["rv", ] <-
      sample(1:3, nr, prob = nodeslst$p[idx_rv], replace = TRUE)

    if ("eta.ka" %in% params) {
      idx_ka <- grep("^eta\\.ka\\.", nodeslst$edge.name)
      rx["eta.ka", ] <-
        sample(0:1, nr, prob = nodeslst$p[idx_ka], replace = TRUE)
    }

    idx_vp <- grep("^eta\\.vp\\.", nodeslst$edge.name)
    idx_q <- grep("^eta\\.q\\.", nodeslst$edge.name)
    idx_vp2 <- grep("^eta\\.vp2\\.", nodeslst$edge.name)
    idx_q2 <- grep("^eta\\.q2\\.", nodeslst$edge.name)
    idx_km <- grep("^eta\\.km\\.", nodeslst$edge.name)

    for (j in 1:nr) {
      cmpt_val <- rx["no.cmpt", j]
      mm_val <- rx["mm", j]

      if (cmpt_val >= 2) {
        rx["eta.vp", j] <- sample(0:1, 1, prob = nodeslst$p[idx_vp])
        rx["eta.q", j] <- sample(0:1, 1, prob = nodeslst$p[idx_q])
      } else {
        rx["eta.vp", j] <- -1
        rx["eta.q", j] <- -1
      }

      if (cmpt_val == 3) {
        rx["eta.vp2", j] <- sample(0:1, 1, prob = nodeslst$p[idx_vp2])
        rx["eta.q2", j] <- sample(0:1, 1, prob = nodeslst$p[idx_q2])
      } else {
        rx["eta.vp2", j] <- -1
        rx["eta.q2", j] <- -1
      }

      if (mm_val == 1) {
        rx["eta.km", j] <- sample(0:1, 1, prob = nodeslst$p[idx_km])
      } else {
        rx["eta.km", j] <- -1
      }
    }

    if (nf > 0) {
      colnames(rx) <- paste0("ant", (nf + 1):nants)
    } else {
      colnames(rx) <- paste0("ant", 1:nants)
    }
  } else {
    rx <- NULL
  }

  if (!is.null(fx) && !is.null(rx)) {
    result <- cbind(fx, rx)
  } else if (!is.null(fx)) {
    result <- fx
  } else {
    result <- rx
  }

  result <- result[params, , drop = FALSE]
  return(result)
}


#' Update pheromone levels for each decision node
#'
#' Compute pheromone increments (\code{delta_phi}) for each node in the ant colony
#' optimization search tree and update the global pheromone levels (\code{phi})
#' based on the ants' paths in the current round.
#'
#' @param r Integer. Current optimization round.
#' @param search.space Character, one of "ivbase" or "oralbase".
#'   Default is "ivbase".
#' @param Q A positive numeric value. Pheromone scaling constant controlling the
#'   amount of pheromone deposited by high-quality solutions during each
#'   iteration. Defaults to 1.
#' @param fitness_history Data frame. History of ants' fitness values and decision
#'   variable selections across rounds.
#' @param nodeslst.hist Data frame. History of node-level pheromone values
#'   across previous rounds.
#' @param alpha A non-negative numeric value. Exponent controlling the influence
#'   of pheromone values on the probability of selecting a component during
#'   solution construction. Defaults to 1.
#' @param rho Numeric in (0, 1). Pheromone evaporation rate. Higher values
#'   increase evaporation, encouraging exploration. Defaults to 0.5.
#' @param diff_tol Numeric. Tolerance threshold controlling when differences
#'   in fitness values are treated as meaningful during pheromone updates.
#'   Defaults to 1.
#' @param phi0 A non-negative numeric value. Initial pheromone value assigned
#' to all nodes at the start of the search. Defaults to 2.
#' @param phi_min A non-negative numeric value. Lower bound for pheromone values, preventing
#'   premature convergence. Defaults to 1.
#' @param phi_max A non-negative numeric value. Upper bound for pheromone values, limiting
#'   excessive reinforcement. Defaults to Inf.
#'
#' @details
#' The update proceeds as follows:
#' \itemize{
#'   \item Initialize the node list for the given search space with \eqn{phi = 0}.
#'   \item Subset the ants from the current round in `fitness_history`.
#'   \item Compute rank-based weights so better-performing ants contribute more:
#'     \deqn{\Delta \phi \propto 1 / \mathrm{rank}^{\alpha}.}
#'   \item Extract the decision columns and attach the computed weights to form
#'     a working table of ant paths and contributions.
#'   \item Map local decision indices to global node numbers using `node.no`
#'     and `local.edge.no` from the node list.
#'   \item For each node, sum contributions from ants that selected the node to
#'     obtain \eqn{\Delta \phi}, then update pheromone with evaporation:
#'     \deqn{\phi_{\mathrm{new}} = (1 - \rho)\,\phi_{\mathrm{prev}} + \Delta \phi.}
#'   \item Clamp updated \eqn{\phi} to be between \code{phi_min} and \code{phi_max}.
#' }
#'
#' @return A data frame (node list) with updated \code{phi} and \code{delta_phi}
#'   for each node.
#'
#' @seealso \link{initNodeList}, \link{rank_new}
#'
#' @author Zhonghui Huang
#'
#' @examples
#' # Define search space
#' search.space <- "ivbase"
#' # Example fitness_history from round 1
#' fitness_history <- data.frame(
#'   round   = rep(1, 8),
#'   mod.no  = 1:8,
#'   no.cmpt = c(1, 1, 2, 2, 3, 3, 2, 2),
#'   eta.km  = c(0, 0, 0, 0, 0, 0, 0, 0),
#'   eta.vc  = c(0, 0, 0, 0, 0, 0, 1, 1),
#'   eta.vp  = c(0, 0, 0, 0, 0, 0, 0, 1),
#'   eta.vp2 = c(0, 0, 0, 0, 0, 0, 0, 0),
#'   eta.q   = c(0, 0, 0, 0, 0, 0, 0, 0),
#'   eta.q2  = c(0, 0, 0, 0, 0, 0, 0, 0),
#'   mm      = c(0, 0, 0, 0, 0, 0, 1, 0),
#'   mcorr   = c(0, 0, 0, 0, 0, 0, 0, 0),
#'   rv      = c(1, 2, 1, 2, 1, 2, 1, 1),
#'   fitness = c(1243.874, 1200.762, 31249.876, 31202.200,
#'               51259.286, 51204.839, 61032.572, 41031.825),
#'   allrank = c(2, 1, 4, 3, 7, 6, 8, 5)
#' )
#'
#' # Example node list history
#' nodeslst.hist <- initNodeList(
#'   search.space = search.space,
#'   phi0 = 2
#' )
#'
#'  phi.calculate(
#'   r = 1,
#'   search.space = search.space,
#'   fitness_history = fitness_history,
#'   nodeslst.hist = nodeslst.hist
#' )
#'
#' @export
phi.calculate <- function(r,
                          search.space = "ivbase",
                          fitness_history = NULL,
                          nodeslst.hist = NULL,
                          Q  = 1,
                          alpha = 1,
                          rho = 0.5,
                          diff_tol = 1,
                          phi0 = 2,
                          phi_min = 1,
                          phi_max = Inf) {

  search.space <-
    match.arg(search.space, choices = c("ivbase", "oralbase", "custom"))
  if (identical(search.space, "custom")) {
    stop(
      "aco currently does not support search.space = 'custom'. Use 'ivbase' or 'oralbase'.",
      call. = FALSE
    )
  }
  node.list <-
    initNodeList(search.space = search.space, phi0 = 0)
  node.list$travel <- r

  current_round_ants <- subset(fitness_history, round == r)

  fitness_history$allrank <-
    rank_new(fitness_history$fitness, diff_tol = diff_tol)

  row.start <- nrow(fitness_history) - nrow(current_round_ants) + 1
  row.stop  <- nrow(fitness_history)

  phi_values <-
    (Q / fitness_history[row.start:row.stop,]$allrank) ^ alpha
  phi_values <- pmax(phi_values, 0)

  decision_cols <-
    grep("no.cmpt|eta|mm|mcorr|rv", names(current_round_ants))
  phi.dat <- as.data.frame(current_round_ants[, decision_cols])
  phi.dat$phi <- phi_values

  if (search.space == "ivbase") {
    col_to_group <- c(
      "no.cmpt" = 1,
      "eta.vp2" = 2,
      "eta.q2"  = 3,
      "eta.vp"  = 4,
      "eta.q"   = 5,
      "eta.vc"  = 6,
      "mm"      = 7,
      "eta.km"  = 8,
      "mcorr"   = 9,
      "rv"      = 10
    )
  } else if (search.space == "oralbase") {
    col_to_group <- c(
      "no.cmpt" = 1,
      "eta.vp2" = 2,
      "eta.q2"  = 3,
      "eta.vp"  = 4,
      "eta.q"   = 5,
      "eta.vc"  = 6,
      "eta.ka"  = 7,
      "mm"      = 8,
      "eta.km"  = 9,
      "mcorr"   = 10,
      "rv"      = 11
    )
  } else {
    stop("Unknown search.space type: must be 'ivbase' or 'oralbase'")
  }

  for (col in setdiff(names(phi.dat), "phi")) {
    group_id <- col_to_group[[col]]
    mapping <-
      stats::setNames(node.list$edge.no[node.list$node.no == group_id],
               node.list$local.edge.no[node.list$node.no == group_id])
    phi.dat[[col]] <- mapping[as.character(phi.dat[[col]])]
  }

  for (n in 1:nrow(node.list)) {
    group_id <- node.list$node.no[n]

    col_name <- names(col_to_group)[col_to_group == group_id]
    if (length(col_name) == 0)
      next

    chosen_mask <- phi.dat[[col_name]] == n
    node.list[n,]$delta_phi <- sum(phi.dat$phi[chosen_mask])
    node.list[n, ]$phi <-
      (1 - rho) * sum(nodeslst.hist[nodeslst.hist$edge.no == n &
                                          nodeslst.hist$travel > (r - 2), ]$phi) + node.list[n, ]$delta_phi

    node.list[n,]$phi <- pmax(node.list[n,]$phi, phi_min)
    node.list[n,]$phi <- pmin(node.list[n,]$phi, phi_max)
  }

  return(node.list)
}


#' Calculate selection probabilities for each node
#'
#' Calculates the probability of selecting each node in an ant colony
#' optimization search, based on pheromone levels \eqn{\phi}.
#'
#' @param nodeslst A data frame of nodes, including columns:
#'   \describe{
#'     \item{phi}{Current pheromone level \eqn{\phi}}
#'     \item{node.no}{Group ID for the decision step}
#'     \item{p}{Probability of selection (to be calculated)}
#'   }
#' @param prob_min Numeric scalar. Minimum probability each node is allowed to have
#'   within its decision group. Set to NULL or 0 to disable smoothing.
#'
#' @details
#' Within each decision group \eqn{G}, selection probabilities are computed from
#' pheromone levels \eqn{\phi} as:
#' \deqn{p_i = \frac{\phi_i}{\sum_{j \in G} \phi_j}}
#'
#' If prob_min is enabled and any calculated probability falls below this value,
#' the algorithm:
#' \enumerate{
#'   \item Sets all probabilities below prob_min to prob_min.
#'   \item Redistributes the remaining probability mass proportionally among
#'         the other nodes in the same group.
#' }
#'
#' This acts as a probability smoothing mechanism, preventing premature
#' convergence by ensuring all nodes retain some chance of being explored.
#'
#' @return The updated node list with recalculated p values.
#'
#' @author Zhonghui Huang
#'
#' @examples
#' node.list <- initNodeList(search.space = "ivbase", phi0 = 1)
#' p.calculation(nodeslst = node.list, prob_min = 0.2)
#' @export
#'
p.calculation <- function(nodeslst,
                          prob_min = NULL) {
  for (group_id in unique(nodeslst$node.no)) {
    group_idx <- which(nodeslst$node.no == group_id)
    group_phi <- nodeslst$phi[group_idx]
    # Base probability calculation
    group_p <- group_phi / sum(group_phi)
    # Apply probability floor if requested
    if (!is.null(prob_min) && prob_min > 0) {
      below_floor <- group_p < prob_min
      if (any(below_floor)) {
        # Fix those below the floor
        group_p[below_floor] <- prob_min
        # After adjustment, check if sum > 1
        total <- sum(group_p)
        if (total > 1) {
          # Too large → scale down the others
          above_floor <- !below_floor
          excess <- total - 1
          group_p[above_floor] <- group_p[above_floor] -
            excess * (group_p[above_floor] / sum(group_p[above_floor]))
        }
      }
    }
    # Store back
    nodeslst$p[group_idx] <- group_p
  }
  return(nodeslst)
}
