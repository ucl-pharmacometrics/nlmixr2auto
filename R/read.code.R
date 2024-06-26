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

read.code2 <- function(search.space,
                      sel.best.code) {
  if (search.space == 1) {
    cmpt.iv = sel.best.code[1]
    eta.km = sel.best.code[2]
    eta.vc = sel.best.code[3]
    eta.vp = sel.best.code[4]
    eta.vp2 = sel.best.code[5]
    eta.q = sel.best.code[6]
    eta.q2 = sel.best.code[7]
    mm = sel.best.code[8]
    mcorr = sel.best.code[9]
    rv = sel.best.code[10]

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
    
    if (cmpt.iv == 1) {
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
    
    
    if (cmpt.iv == 2) {
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
    
    
    if (cmpt.iv == 3) {
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
    if (rv == 1) {
      mod.lib10 <- ",additive"
    }

    if (rv == 2) {
      mod.lib10 <- ",proportional"
    }
    if (rv == 3) {
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

















