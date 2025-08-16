#' Create Ant Population for Ant Colony Optimization (ACO)
#'
#' This function generates a population of "ants" (candidate model encodings)
#' for use in the ACO algorithm. Each ant is represented as a vector of binary
#' or categorical choices corresponding to structural model components
#' (e.g., compartments, random effects, error models).
#'
#' Ants can be generated either by:
#' \itemize{
#'   \item \strong{Fixed initialization} — the first 6 ants correspond to simple,
#'         pre-defined base models (1–3 compartments, with/without variance terms).
#'   \item \strong{Probability sampling} — additional ants are sampled from the
#'         pheromone-based probability distributions stored in \code{node.list}.
#' }
#'
#' @param search.space Character string. Defines the search space type.
#'   Options are:
#'   \itemize{
#'     \item `"ivbase"` — IV models (default).
#'     \item `"oralbase"` — Oral models (includes \code{eta.ka} as a parameter).
#'   }
#' @param no.ants Integer. Number of ants (candidate models) to generate.
#'   Must be at least 6.
#' @param initialize Logical. If \code{TRUE}, the first 6 ants are fixed to
#'   simple base models, and any remaining ants are sampled probabilistically.
#'   If \code{FALSE}, all ants are sampled probabilistically.
#' @param node.list A data structure (usually from \code{initNodeList()})
#'   containing node names and their associated sampling probabilities
#'   (pheromone levels). Used to guide random sampling of ants.
#'
#' @return A numeric matrix where:
#' \itemize{
#'   \item Rows correspond to structural/model features (\code{"no.cmpt"},
#'         \code{"eta.km"}, \code{"eta.vc"}, \code{"eta.vp"}, \code{"eta.vp2"},
#'         \code{"eta.q"}, \code{"eta.q2"}, optionally \code{"eta.ka"},
#'         \code{"mm"}, \code{"mcorr"}, \code{"rv"}).
#'   \item Columns correspond to individual ants (\code{"ant1"}, \code{"ant2"}, ...).
#' }
#' Values are integers (0/1 for inclusion/exclusion, -1 for "not applicable",
#' categorical codes for model features).
#'
#' @examples
#' \dontrun{
#' # Example: initialize 10 ants in IV search space
#' node.list <- initNodeList(search.space = "ivbase", initial.phi = 0.5)
#' ants <- createAnts(search.space = "ivbase", no.ants = 10,
#'                    initialize = TRUE, node.list = node.list)
#' ants
#' }
#'
#' @seealso \code{\link{aco.operator}}, \code{\link{initNodeList}}
#'
#' @export

createAnts <- function(search.space="ivbase",
                       no.ants=6,
                       initialize=F,
                       node.list=NULL) {
  # --- Defensive check: no.ants >= 6 ---
  if (no.ants < 6) {
    stop("Number of ants (no.ants) must be at least 6.")
  }

  # --- Determine if eta.ka is present (oralbase) ---
  has_eta_ka <- any(grepl("^eta\\.ka", node.list$node.names))

  # --- Helper: get index in node.list by name ---
  idx <- function(names)
    which(node.list$node.names %in% names)

  # --- Prepare empty storage ---
  no.cmpt <- data.frame(matrix(nrow = 1, ncol = no.ants))
  eta.km  <- data.frame(matrix(nrow = 1, ncol = no.ants))
  eta.vc  <- data.frame(matrix(nrow = 1, ncol = no.ants))
  eta.vp  <- data.frame(matrix(nrow = 1, ncol = no.ants))
  eta.vp2 <- data.frame(matrix(nrow = 1, ncol = no.ants))
  eta.q   <- data.frame(matrix(nrow = 1, ncol = no.ants))
  eta.q2  <- data.frame(matrix(nrow = 1, ncol = no.ants))
  if (has_eta_ka) {
    eta.ka <- data.frame(matrix(nrow = 1, ncol = no.ants))
  }
  mm      <- data.frame(matrix(nrow = 1, ncol = no.ants))
  mcorr   <- data.frame(matrix(nrow = 1, ncol = no.ants))
  rv      <- data.frame(matrix(nrow = 1, ncol = no.ants))

  # --- Initialize population ---
  if (initialize) {
    # First 6 ants: fixed simple model
    no.cmpt[1, 1:6]  <- c(1, 1, 2, 2, 3, 3)
    eta.vc[1, 1:6]   <- 0
    eta.vp[1, 1:6]   <- c(-1, -1, 0, 0, 0, 0)
    eta.q[1, 1:6]    <- c(-1, -1, 0, 0, 0, 0)
    eta.vp2[1, 1:6]  <- c(-1, -1, -1, -1, 0, 0)
    eta.q2[1, 1:6]   <- c(-1, -1, -1, -1, 0, 0)
    rv[1, 1:6]       <- c(1, 2, 1, 2, 1, 2)
    mcorr[1, 1:6]    <- 0
    mm[1, 1:6]       <- 0
    eta.km[1, 1:6]   <- -1
    if (has_eta_ka)
      eta.ka[1, 1:6] <- -1

    # Remaining ants: probability sampling
    if (no.ants > 6) {
      no.cmpt[1, 7:no.ants] <-
        sample(1:3, no.ants - 6, prob = node.list[idx(c("1Cmpt", "2Cmpt", "3Cmpt")),]$p, replace = TRUE)
      eta.vc[1, 7:no.ants]  <-
        sample(0:1, no.ants - 6, prob = node.list[idx(c("eta.vc.no", "eta.vc.yes")),]$p, replace = TRUE)
      mm[1, 7:no.ants]      <-
        sample(0:1, no.ants - 6, prob = node.list[idx(c("mm.no", "mm.yes")),]$p, replace = TRUE)
      mcorr[1, 7:no.ants]   <-
        sample(0:1, no.ants - 6, prob = node.list[idx(c("mcorr.no", "mcorr.yes")),]$p, replace = TRUE)
      rv[1, 7:no.ants]      <-
        sample(1:3, no.ants - 6, prob = node.list[idx(c("add", "prop", "comb")),]$p, replace = TRUE)
      if (has_eta_ka) {
        eta.ka[1, 7:no.ants] <-
          sample(0:1,
                 no.ants - 6,
                 prob = node.list[idx(c("eta.ka.no", "eta.ka.yes")),]$p,
                 replace = TRUE)
      }

      for (j in 7:no.ants) {
        if (no.cmpt[1, j] == 3) {
          eta.vp2[1, j] <-
            sample(0:1, 1, prob = node.list[idx(c("eta.vp2.no", "eta.vp2.yes")),]$p)
          eta.q2[1, j]  <-
            sample(0:1, 1, prob = node.list[idx(c("eta.q2.no", "eta.q2.yes")),]$p)
        } else {
          eta.vp2[1, j] <- -1
          eta.q2[1, j]  <- -1
        }

        if (no.cmpt[1, j] > 1) {
          eta.vp[1, j] <-
            sample(0:1, 1, prob = node.list[idx(c("eta.vp.no", "eta.vp.yes")),]$p)
          eta.q[1, j]  <-
            sample(0:1, 1, prob = node.list[idx(c("eta.q.no", "eta.q.yes")),]$p)
        } else {
          eta.vp[1, j] <- -1
          eta.q[1, j]  <- -1
        }

        if (mm[1, j] == 1) {
          eta.km[1, j] <-
            sample(0:1, 1, prob = node.list[idx(c("eta.km.no", "eta.km.yes")),]$p)
        } else {
          eta.km[1, j] <- -1
        }
      }
    }
  } else {
    # Fully sampled initialization
    no.cmpt[1, ] <-
      sample(1:3, no.ants, prob = node.list[idx(c("1Cmpt", "2Cmpt", "3Cmpt")),]$p, replace = TRUE)
    eta.vc[1, ]  <-
      sample(0:1, no.ants, prob = node.list[idx(c("eta.vc.no", "eta.vc.yes")),]$p, replace = TRUE)
    mm[1, ]      <-
      sample(0:1, no.ants, prob = node.list[idx(c("mm.no", "mm.yes")),]$p, replace = TRUE)
    mcorr[1, ]   <-
      sample(0:1, no.ants, prob = node.list[idx(c("mcorr.no", "mcorr.yes")),]$p, replace = TRUE)
    rv[1, ]      <-
      sample(1:3, no.ants, prob = node.list[idx(c("add", "prop", "comb")),]$p, replace = TRUE)
    if (has_eta_ka) {
      eta.ka[1, ] <-
        sample(0:1, no.ants, prob = node.list[idx(c("eta.ka.no", "eta.ka.yes")),]$p, replace = TRUE)
    }

    for (j in 1:no.ants) {
      if (no.cmpt[1, j] == 3) {
        eta.vp2[1, j] <-
          sample(0:1, 1, prob = node.list[idx(c("eta.vp2.no", "eta.vp2.yes")),]$p)
        eta.q2[1, j]  <-
          sample(0:1, 1, prob = node.list[idx(c("eta.q2.no", "eta.q2.yes")),]$p)
      } else {
        eta.vp2[1, j] <- -1
        eta.q2[1, j]  <- -1
      }

      if (no.cmpt[1, j] > 1) {
        eta.vp[1, j] <-
          sample(0:1, 1, prob = node.list[idx(c("eta.vp.no", "eta.vp.yes")),]$p)
        eta.q[1, j]  <-
          sample(0:1, 1, prob = node.list[idx(c("eta.q.no", "eta.q.yes")),]$p)
      } else {
        eta.vp[1, j] <- -1
        eta.q[1, j]  <- -1
      }

      if (mm[1, j] == 1) {
        eta.km[1, j] <-
          sample(0:1, 1, prob = node.list[idx(c("eta.km.no", "eta.km.yes")),]$p)
      } else {
        eta.km[1, j] <- -1
      }
    }
  }

  # --- Combine into final matrix ---
  if (has_eta_ka) {
    ants.all <-
      rbind(no.cmpt,
            eta.km,
            eta.vc,
            eta.vp,
            eta.vp2,
            eta.q,
            eta.q2,
            eta.ka,
            mm,
            mcorr,
            rv)
    ants.all <- as.matrix(ants.all)
    colnames(ants.all) <- paste0("ant", seq_len(ncol(ants.all)))
    rownames(ants.all) <- c(
      "no.cmpt", "eta.km", "eta.vc", "eta.vp", "eta.vp2",
      "eta.q", "eta.q2", "eta.ka", "mm", "mcorr", "rv"
    )
  } else {
    ants.all <-
      rbind(no.cmpt,
            eta.km,
            eta.vc,
            eta.vp,
            eta.vp2,
            eta.q,
            eta.q2,
            mm,
            mcorr,
            rv)
    ants.all <- as.matrix(ants.all)
    colnames(ants.all) <- paste0("ant", seq_len(ncol(ants.all)))
    rownames(ants.all) <- c(
      "no.cmpt", "eta.km", "eta.vc", "eta.vp", "eta.vp2",
      "eta.q", "eta.q2", "mm", "mcorr", "rv"
    )
  }
  return(ants.all)
}
