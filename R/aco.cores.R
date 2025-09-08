#' Update pheromone levels for each decision node
#'
#' Computes the pheromone increment (\eqn{\Delta \phi}) for each node
#' in the ant colony optimization (ACO) search tree, and updates the global pheromone
#' levels (\eqn{\phi}) based on the ants' paths in a given round.
#'
#' @param r Integer. The current optimization round.
#' @param search.space Character. Either `"ivbase"` or `"oralbase"`, determines
#'   which decision variables and node groups are used.
#' @param fitness_history Data frame. History of all ants' fitness values and decision
#'   variable selections across rounds.
#' @param node.list.history Data frame. History of node-level pheromone values
#'   across previous rounds.
#' @param alpha.value Numeric. Exponent used when converting ranks to pheromone
#'   contributions: \eqn{\Delta\phi \propto 1/\text{rank}^\alpha}.
#' @param rho Numeric. Pheromone evaporation rate (fraction evaporated per iteration).
#' @param sig.diff Numeric. Significance threshold for distinguishing fitness ranks.
#' @param lower.limit.phi Numeric. Lower bound on updated \eqn{\phi}.
#' @param upper.limit.phi Numeric. Upper bound on updated \eqn{\phi}.
#'
#' @details
#' Steps performed:
#' 1. **Initialize node list** for the given search space with \eqn{\phi = 0}.
#' 2. **Subset ants of current round** from `fitness_history`.
#' 3. **Compute rank-based \eqn{\Delta\phi} weights**: ants with higher fitness
#'    receive larger contributions, scaled by `alpha.value`.
#' 4. **Extract decision columns** (variables controlling path choices) and append
#'    the computed \eqn{\phi} values as a new column in `phi.dat`.
#' 5. **Map local decision indices to global node numbers** using `node.group` and
#'    `local.node.no` from `node.list`.
#'    - For example, in `"ivbase"`, `"eta.km"` maps to group 8; in `"oralbase"`, `"eta.ka"`
#'      is included as group 7.
#' 6. **For each node**:
#'    - Identify ants that chose this node (`chosen_mask`).
#'    - Sum their \eqn{\phi} values to get `delta_phi`.
#'    - Update global \eqn{\phi} with evaporation:
#'      \deqn{\phi_{\text{new}} = (1-\rho) \cdot \phi_{\text{prev, recent}} + \Delta\phi}
#'    - Clip \eqn{\phi} to `[lower.limit.phi, upper.limit.phi]`.
#'
#' This update rule maintains both exploration (through evaporation) and exploitation
#' (reinforcing paths chosen by better-performing ants).
#'
#' @return A data frame (node list) with updated `phi` and `delta_phi` for each node.
#'
#' @seealso \code{\link{initNodeList}}, \code{\link{rank_new}}
#'
#' @examples
#' # Define search space
#' search.space <- "ivbase"
#' initial.phi <- 0
#'
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
#' node.list.history <- initNodeList(
#'   search.space = search.space,
#'   initial.phi = initial.phi
#' )
#'
#' # Call phi.calculate
#' updated_nodes <- phi.calculate(
#'   r = 1,
#'   search.space = search.space,
#'   fitness_history = fitness_history,
#'   node.list.history = node.list.history,
#'   alpha.value = 1,
#'   rho = 0.2,
#'   sig.diff = 1,
#'   lower.limit.phi = 1,
#'   upper.limit.phi = Inf
#' )
#' print(updated_nodes )
#'
#' @export
phi.calculate <- function(r,
                          search.space = "ivbase",
                          fitness_history = NULL,
                          node.list.history = NULL,
                          param.Q  = 1,
                          alpha.value = 1,
                          rho = 0.2,
                          sig.diff = 1,
                          lower.limit.phi = 1,
                          upper.limit.phi = Inf) {
  node.list <-
    initNodeList(search.space = search.space, initial.phi = 0)
  node.list$travel <- r

  current_round_ants <- subset(fitness_history, round == r)

  fitness_history$allrank <-
    rank_new(fitness_history$fitness, sig.diff = sig.diff)

  row.start <- nrow(fitness_history) - nrow(current_round_ants) + 1
  row.stop  <- nrow(fitness_history)

  phi_values <-
    (param.Q / fitness_history[row.start:row.stop,]$allrank) ^ alpha.value
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
    mapping <- setNames(node.list$node.no[node.list$node.group == group_id],
                        node.list$local.node.no[node.list$node.group == group_id])
    phi.dat[[col]] <- mapping[as.character(phi.dat[[col]])]
  }

  for (n in 1:nrow(node.list)) {
    group_id <- node.list$node.group[n]

    col_name <- names(col_to_group)[col_to_group == group_id]
    if (length(col_name) == 0)
      next

    chosen_mask <- phi.dat[[col_name]] == n

    node.list[n,]$delta_phi <- sum(phi.dat$phi[chosen_mask])

    node.list[n,]$phi <-
      (1 - rho) * sum(node.list.history[node.list.history$node.no == n &
                                          node.list.history$travel > (r - 2),]$phi) + node.list[n,]$delta_phi

    node.list[n,]$phi <- pmax(node.list[n,]$phi, lower.limit.phi)
    node.list[n,]$phi <- pmin(node.list[n,]$phi, upper.limit.phi)
  }

  return(node.list)
}




#' Calculate selection probabilities for each node
#'
#' This function calculates the probability (`p`) of selecting each node in an
#' Ant Colony Optimization (ACO) search, based on current pheromone levels (`phi`).
#' It supports an optional probability floor (`prob.floor`) mechanism that ensures
#' no node's probability drops below a minimum threshold, redistributing the
#' remaining probability proportionally among other nodes in the same group.
#'
#' @param search.space Character string or numeric code indicating the search space.
#'   Accepts `"ivbase"`, `"oralbase"`, or `1` (equivalent to `"ivbase"`).
#' @param node.list A data frame of nodes, including columns:
#'   \describe{
#'     \item{phi}{Current pheromone level}
#'     \item{node.group}{Group ID for the decision step}
#'     \item{p}{Probability of selection (to be calculated)}
#'   }
#' @param prob.floor Numeric scalar. Minimum probability each node is allowed to have
#'   within its decision group. Set to `NULL` or `0` to disable smoothing.
#'
#' @details
#' The probability for each node in a group is calculated as:
#' \deqn{p_i = \frac{\phi_i}{\sum_{j \in G} \phi_j}}
#'
#' If \code{prob.floor} is set and any calculated probability falls below this value,
#' the algorithm:
#' \enumerate{
#'   \item Sets all probabilities below \code{prob.floor} to \code{prob.floor}.
#'   \item Redistributes the remaining probability mass proportionally among
#'         the other nodes in the same group.
#' }
#'
#' This acts as a \strong{probability smoothing mechanism}, preventing premature
#' convergence by ensuring all nodes retain some chance of being explored.
#'
#' @return The updated \code{node.list} with recalculated \code{p} values.
#'
#' @examples
#' \dontrun{
#' node.list <- initNodeList(search.space = "ivbase", initial.phi = 1)
#' updated.nodes <- p.calculation(
#'                               node.list =  node.list,
#'                               prob.floor = 0.2)
#' }
#' @export
#'
p.calculation <- function(node.list,
                          prob.floor = NULL) {
  for (group_id in unique(node.list$node.group)) {
    group_idx <- which(node.list$node.group == group_id)
    group_phi <- node.list$phi[group_idx]

    # Base probability calculation
    group_p <- group_phi / sum(group_phi)

    # Apply probability floor if requested
    if (!is.null(prob.floor) && prob.floor > 0) {
      below_floor <- group_p < prob.floor

      if (any(below_floor)) {
        # Fix those below the floor
        group_p[below_floor] <- prob.floor

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
    node.list$p[group_idx] <- group_p
  }
  return(node.list)
}
