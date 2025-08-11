#' Generate pharmacokinetic model name based on selected codes
#'
#' Generate a model name based on the selected best code and specified search space.
#'
#' @param search.space An integer indicating the search space.
#' @param sel.best.code A numeric vector containing the selected best codes. If `search.space` is 1, it should have 12 elements representing different components and parameters.
#'
#' @details The function uses the values in `sel.best.code` to determine the components and parameters of the model name. The model name is constructed based on the following:
#' \itemize{
#'   \item \code{cmpt.iv1} and \code{cmpt.iv2}: Determine the number of compartments (1, 2, or 3).
#'   \item \code{eta.km}, \code{eta.vc}, \code{eta.vp}, \code{eta.vp2}, \code{eta.q}, and \code{eta.q2}: Determine the presence of various random effects.
#'   \item \code{mm}: Determine the elimination method (first-order or Michaelis-Menten).
#'   \item \code{mcorr}: Determine whether correlation is present.
#'   \item \code{rv}: Determine the residual variance type (additive, proportional, or combined).
#' }
#'
#' @return A character string representing the model name.
#'
#' @examples
#' sel.best.code <- c(1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1)
#' model_name <- read.code(1, sel.best.code)
#' print(model_name)
#'
#' @export

read.code <- function(search.space,
                      sel.best.code) {
  if (search.space == 1) {
    cmpt.iv1 = sel.best.code[1]
    cmpt.iv2 = sel.best.code[2]
    eta.km = sel.best.code[3]
    eta.vc = sel.best.code[4]
    eta.vp = sel.best.code[5]
    eta.vp2 = sel.best.code[6]
    eta.q = sel.best.code[7]
    eta.q2 = sel.best.code[8]
    mm = sel.best.code[9]
    mcorr = sel.best.code[10]
    rv1 = sel.best.code[11]
    rv2 = sel.best.code[12]

    mod.lib1 <- NULL
    mod.lib2 <- NULL
    mod.lib3 <- NULL
    mod.lib4 <- NULL
    mod.lib5 <- NULL
    mod.lib6 <- NULL
    mod.lib7 <- NULL
    mod.lib8 <- NULL
    mod.lib9 <- NULL
    mod.lib10 <- NULL

    if (cmpt.iv1 == 0 & cmpt.iv2 == 0) {
      mod.lib1 <- "1Cmpt,IIV.cl"

      if (eta.vc == 0) {
        mod.lib2 <- ""
      }
      if (eta.vc == 1) {
        mod.lib2 <- ".vc"
      }
      if (eta.km == 0) {
        mod.lib7 <- ""
      }
      if (eta.km == 1) {
        mod.lib7 <- ".km"
      }
    }

    if (cmpt.iv1 == 0 & cmpt.iv2 == 1) {
      mod.lib1 <- "1Cmpt,IIV.cl"

      if (eta.vc == 0) {
        mod.lib2 <- ""
      }
      if (eta.vc == 1) {
        mod.lib2 <- ".vc"
      }
      if (eta.km == 0) {
        mod.lib7 <- ""
      }
      if (eta.km == 1) {
        mod.lib7 <- ".km"
      }
    }

    if (cmpt.iv1 == 1 & cmpt.iv2 == 0) {
      mod.lib1 <- "2Cmpt,IIV.cl"

      if (eta.vc == 0) {
        mod.lib2 <- ""
      }
      if (eta.vc == 1) {
        mod.lib2 <- ".vc"
      }

      if (eta.vp == 0) {
        mod.lib3 <- ""
      }
      if (eta.vp == 1) {
        mod.lib3 <- ".vp"
      }

      if (eta.q == 0) {
        mod.lib4 <- ""
      }
      if (eta.q == 1) {
        mod.lib4 <- ".q"
      }
      if (eta.km == 0) {
        mod.lib7 <- ""
      }
      if (eta.km == 1) {
        mod.lib7 <- ".km"
      }
    }


    if (cmpt.iv1 == 1 & cmpt.iv2 == 1) {
      mod.lib1 <- "3Cmpt,IIV.cl"


      if (eta.vc == 0) {
        mod.lib2 <- ""
      }
      if (eta.vc == 1) {
        mod.lib2 <- ".vc"
      }

      if (eta.vp == 0) {
        mod.lib3 <- ""
      }
      if (eta.vp == 1) {
        mod.lib3 <- ".vp"
      }

      if (eta.vp2 == 0) {
        mod.lib4 <- ""
      }
      if (eta.vp2 == 1) {
        mod.lib4 <- ".vp2"
      }
      if (eta.q == 0) {
        mod.lib5 <- ""
      }
      if (eta.q == 1) {
        mod.lib5 <- ".q"
      }
      if (eta.q2 == 0) {
        mod.lib6 <- ""
      }
      if (eta.q2 == 1) {
        mod.lib6 <- ".q2"
      }
      if (eta.km == 0) {
        mod.lib7 <- ""
      }
      if (eta.km == 1) {
        mod.lib7 <- ".km"
      }
    }


    if (mm == 0) {
      mod.lib8 <- ",first-order_elminiation"
    }
    if (mm == 1) {
      mod.lib8 <- ",M-M_elimination"
    }

    if (mcorr == 0) {
      mod.lib9 <- ",nocorr"
    }
    if (mcorr == 1) {
      mod.lib9 <- ",full_omega_matrix"
    }
    if (rv1 == 0 & rv2 == 0) {
      mod.lib10 <- ",additive"
    }
    if (rv1 == 0 & rv2 == 1) {
      mod.lib10 <- ",additive"
    }
    if (rv1 == 1 & rv2 == 0) {
      mod.lib10 <- ",proportional"
    }
    if (rv1 == 1 & rv2 == 1) {
      mod.lib10 <- ",combined"
    }
    mod.name <-
      paste(
        mod.lib1,
        mod.lib2,
        mod.lib3,
        mod.lib4,
        mod.lib5,
        mod.lib6,
        mod.lib7,
        mod.lib8,
        mod.lib9,
        mod.lib10
      )
    cleaned_mod.name <- gsub("\\s+", "",  mod.name)

    return(cleaned_mod.name)

  }
}





#' Generate pharmacokinetic model name based on selected codes [2]
#'
#' Generate a model name based on the selected best code and specified search space for categorical coding algorithm
#'
#' @param search.space An integer indicating the search space.
#' @param sel.best.code A numeric vector containing the selected best codes. If `search.space` is 1, it should have 12 elements representing different components and parameters.
#'
#' @details The function uses the values in `sel.best.code` to determine the components and parameters of the model name. The model name is constructed based on the following:
#' \itemize{
#'   \item \code{cmpt}: Determine the number of compartments (1, 2, or 3).
#'   \item \code{eta.km}, \code{eta.vc}, \code{eta.vp}, \code{eta.vp2}, \code{eta.q}, and \code{eta.q2}: Determine the presence of various random effects.
#'   \item \code{mm}: Determine the elimination method (first-order or Michaelis-Menten).
#'   \item \code{mcorr}: Determine whether correlation is present.
#'   \item \code{rv}: Determine the residual variance type (additive, proportional, or combined).
#' }
#'
#' @return A character string representing the model name.
#'
#' @examples
#' sel.best.code <- c(1, 1, 0, 0, 0, 0, 0, 0, 0, 1)
#' model_name <- read.code(1, sel.best.code)
#' print(model_name)
#'
#' @export

CodetoMod <- function(search.space, sel.best.code) {
  # Configuration for each search space
  cfg <- switch(
    search.space,
    "ivbase" = list(
      prefix       = "iv",
      no_cmpt_idx  = 1L,
      # number of compartments
      eta_idx      = 2:7,
      # KM, VC, VP, VP2, Q, Q2
      eta_labels   = c("KM", "VC", "VP", "VP2", "Q", "Q2"),
      mm_idx       = 8L,
      # mm flag
      mcorr_idx    = 9L,
      # correlation flag
      rv_idx       = 10L               # residual error model
    ),
    "oralbase" = list(
      prefix       = "oral",
      no_cmpt_idx  = 1L,
      eta_idx      = 2:8,
      # KM, VC, VP, VP2, Q, Q2, KA
      eta_labels   = c("KM", "VC", "VP", "VP2", "Q", "Q2", "KA"),
      mm_idx       = 9L,
      mcorr_idx    = 10L,
      rv_idx       = 11L
    ),
    stop("Unknown search.space: must be 'ivbase' or 'oralbase'.")
  )

  # Safe integer accessor
  grab_int <- function(i) {
    x <- suppressWarnings(as.integer(sel.best.code[i]))
    if (is.na(x))
      0L
    else
      x
  }

  # Extract fields
  no_cmpt   <- grab_int(cfg$no_cmpt_idx)
  eta_flags <-
    as.integer(sel.best.code[cfg$eta_idx])
  eta_flags[is.na(eta_flags)] <- 0L
  mm        <- grab_int(cfg$mm_idx)
  mcorr     <- grab_int(cfg$mcorr_idx)
  rv        <- grab_int(cfg$rv_idx)

  # Structural string
  cmpt_string <- paste0(no_cmpt, "cmpt")

  # Forced ETA label by mm, plus any active ETA labels
  forced <-
    if (mm == 1L)
      "Vmax"
  else
    "CL"   # mm=1 -> Vmax, mm=0 -> CL
  active_labels <- cfg$eta_labels[which(eta_flags == 1L)]
  eta_parts <- unique(c(forced, active_labels))
  eta_string <- paste0("eta", paste(eta_parts, collapse = ""))

  mm_string <- if (mm == 0L) {
    "First-order elimination"
  } else {
    "Michaelis-Menten elimination"
  }

  corr_string <- if (mcorr == 0L) {
    "No correlation"
  } else {
    "Eta_correlated"
  }

  rv_string <- switch(
    as.character(rv),
    "1" = "additive",
    "2" = "proportional",
    "3" = "combined",
    paste0("rv", rv) # fallback
  )


  # Final model name
  paste(cfg$prefix,
        cmpt_string,
        eta_string,
        mm_string,
        corr_string,
        rv_string,
        sep = "_")
}
