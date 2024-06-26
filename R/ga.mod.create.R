#' Create chromosome string of parameters
#'
#' Converts binary codes into categorical codes and assembles them into a string
#' of model parameters for use in a genetic algorithm. It applies specific rules for the conversion
#' of input variables and calls the external model creation function `ex.mod.create.iv.mm`.
#'
#' @param modi The current model number.
#' @param cmpt.iv1 A binary value representing the first part of the number of compartments.
#' @param cmpt.iv2 A binary value representing the second part of the number of compartments.
#' \itemize{
#'   \item \code{cmpt.iv1 = 0, cmpt.iv2 = 0}: One compartment
#'   \item \code{cmpt.iv1 = 0, cmpt.iv2 = 1}: One compartment
#'   \item \code{cmpt.iv1 = 1, cmpt.iv2 = 0}: Two compartments
#'   \item \code{cmpt.iv1 = 1, cmpt.iv2 = 1}: Three compartments
#' }
#' @param eta.km A binary value indicating the presence (1) or absence (0) of random effects for the Michaelis-Menten constant (Km).
#' @param eta.vc A binary value indicating the presence (1) or absence (0) of random effects for the central volume of distribution (Vc).
#' @param eta.vp A binary value indicating the presence (1) or absence (0) of random effects for the peripheral volume of distribution (Vp).
#' @param eta.vp2 A binary value indicating the presence (1) or absence (0) of random effects for the second peripheral volume of distribution (Vp2).
#' @param eta.q A binary value indicating the presence (1) or absence (0) of random effects for the intercompartmental clearance (Q).
#' @param eta.q2 A binary value indicating the presence (1) or absence (0) of random effects for the second intercompartmental clearance (Q2).
#' @param mm A binary value indicating the presence (1) or absence (0) of Michaelis-Menten elimination.
#' @param mcorr A binary value indicating the presence (1) or absence (0) of correlation between subjects
#' @param rv1 A binary value representing the first part of the residual error.
#' @param rv2 A binary value representing the second part of the residual error.
#'  \itemize{
#'   \item \code{rv1 = 0, rv2 = 0}: Additive residual error
#'   \item \code{rv1 = 0, rv2 = 1}: Additive residual error
#'   \item \code{rv1 = 1, rv2 = 0}: Proportional residual error
#'   \item \code{rv1 = 1, rv2 = 1}: Combined residual error
#' }
#' @return A vector containing the assembled string of model parameters.
#' @import stringr
#'
#' @examples
#' # Example usage:
#' ga.mod.create(modi=1001,cmpt.iv1 = 1, cmpt.iv2 = 1, eta.km = 0, eta.vc = 0,
#'                             eta.vp = 0, eta.vp2 = 1, eta.q = 0, eta.q2 = 0,
#'                             mm = 1, mcorr = 1, rv1 = 1, rv2 = 0)
#'
#' @export

ga.mod.create <- function(modi,
                          cmpt.iv1,
                          cmpt.iv2,
                          eta.km,
                          eta.vc,
                          eta.vp,
                          eta.vp2,
                          eta.q,
                          eta.q2,
                          mm,
                          mcorr,
                          rv1,
                          rv2) {
  # if un.code function is used then, there is no the first situation
  if (cmpt.iv1 == 0 & cmpt.iv2 == 0) {
    cmpt.iv = 1
  }

  if (cmpt.iv1 == 0 & cmpt.iv2 == 1) {
    cmpt.iv = 1
  }

  if (cmpt.iv1 == 1 & cmpt.iv2 == 0) {
    cmpt.iv = 2
  }

  if (cmpt.iv1 == 1 & cmpt.iv2 == 1) {
    cmpt.iv = 3
  }

  # if un.code function is used then, there is no the first situation

  if (rv1 == 0 & rv2 == 0) {
    rv = 1
  }

  if (rv1 == 0 & rv2 == 1) {
    rv = 1
  }

  if (rv1 == 1 & rv2 == 0) {
    rv = 2
  }

  if (rv1 == 1 & rv2 == 1) {
    rv = 3
  }

  ex.mod.create.iv.mm(
    modi = modi,
    cmpt.iv = cmpt.iv,
    eta.km = eta.km,
    eta.vc = eta.vc,
    eta.vp = eta.vc,
    eta.vp2 = eta.vp2,
    eta.q = eta.q,
    eta.q2 = eta.q2,
    mm = mm,
    mcorr = mcorr,
    rv = rv
  )

  string <- c(cmpt.iv,
              eta.km,
              eta.vc,
              eta.vp,
              eta.vp2,
              eta.q,
              eta.q2,
              mm,
              mcorr,
              rv)

  return(string)
}
