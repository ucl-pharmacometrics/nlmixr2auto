#' Generate neighbors data frame
#'
#' Generate a data frame containing all possible neighbor models of a given current model string.
#'
#' @param current_string A named vector or list representing the current state, 
#' where names correspond to variable names and values correspond to their current values.
#'
#' @return A DataFrame containing all possible neighboring states.
#'
#' @examples
#' current_string <- list(cmpt.iv = 2, eta.km = 0, eta.vc = 1, eta.vp = 0, 
#'                       eta.vp2 = 1, eta.q = 0, eta.q2 = 1, mm = 0, 
#'                       mcorr = 1, rv = 2)
#' neighbors_df <- generate_neighbors_df(current_string)
#' print(neighbors_df)
#'
#' @export


generate_neighbors_df <- function(current_string) {
  
  options_list <- list(
    cmpt.iv = 1:3,
    eta.km = 0:1,
    eta.vc = 0:1,
    eta.vp = 0:1,
    eta.vp2 = 0:1,
    eta.q = 0:1,
    eta.q2 = 0:1,
    mm= 0:1,
    mcorr = 0:1,
    rv = 1:3
  )
  
  neighbors <- data.frame(matrix(ncol = length(current_string), nrow = 0))
  colnames(neighbors) <- names(current_string)
  
  for (variable in names(current_string)) {
    options <- options_list[[variable]]
    for (option in options) {
      if (option != current_string[variable]) {
        new_neighbor <- current_string
        new_neighbor[variable] <- option
        neighbors <- rbind(neighbors, new_neighbor)
      }
    }
  }
  return(neighbors)
}