#' Generate initial population for genetic algorithm
#'
#' Generates the initial population of chromosomes for a genetic algorithm.
#' Each chromosome is represented as a binary string with bits that are either 0 or 1.
#'
#' @param npopsize An integer representing the number of individuals (chromosomes) in the population.
#' @param nbits An integer representing the number of bits in each chromosome.
#'
#' @return A matrix of size \code{npopsize} by \code{nbits}, where each element is either 0 or 1.
#'
#' @examples
#' # Create an initial population of 10 chromosomes, each with 12 bits
#' create.pop(npopsize = 10, nbits = 12)
#'
#' @export

create.pop <- function(npopsize,
                       nbits) {
  population <- matrix(as.double(NA),
                       nrow = npopsize,
                       ncol = nbits)
  
  for (j in 1:npopsize) {
    population[j,] <- sample(
      x = c(0, 1),
      size = nbits,
      prob = c(0.5, 0.5),
      replace = T
    )
  }
  
  return(population)
  
}
