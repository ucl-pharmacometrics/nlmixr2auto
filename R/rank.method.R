#' Ranking with significance difference threshold
#'
#' Performs a custom ranking of a numeric vector,and ajusts the ranks of values that
#' differ by less than a specified threshold, ensuring they receive the same rank.
#'
#' @param x1 A numeric vector to be ranked.
#' @param diff_tol A numeric value specifying the significance difference threshold.
#' Values within this threshold are considered equal and receive the same rank.
#'
#' @return A numeric vector representing the adjusted ranks of the input values.
#'
#' @author Zhonghui Huang
#'
#' @examples
#' x1 <- c(10, 20, 20.5, 30)
#' diff_tol <- 1
#' ranked_list <- rank_new(x1, diff_tol)
#' print(ranked_list)
#'
#' @export

rank_new <- function(x1, diff_tol)
{
  rx1 <- rank(x1, ties.method = "min")
  datarank <- data.frame(x1, rx1)
  datarank$id <- seq(1, length(x1), 1)
  datarank2 <- datarank[order(datarank$rx1), ]
  datarank2$rx2 <- NA
  datarank2[1, ]$rx2 <- 1

  for (i in 2:nrow(datarank2)) {
    test = datarank2[i, ]$x1 - datarank2[i - 1, ]$x1
    if (isTRUE(test < diff_tol) == TRUE) {
      datarank2[i, ]$rx2 = datarank2[i - 1, ]$rx2
    }
    else
    {
      datarank2[i, ]$rx2 = datarank2[i, ]$rx1
    }
  }

  datarank3 <- datarank2[order(datarank2$id), ]
  ranklist <- datarank3$rx2
  return(ranklist)

}
