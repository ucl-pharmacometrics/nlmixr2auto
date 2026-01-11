#' Create an initial GA population
#'
#' Generates an initial population for a genetic algorithm (GA). Each individual
#' is a binary chromosome represented by a numeric vector containing 0 and 1.
#'
#' @param npop Integer. Number of individuals (chromosomes) in the population.
#' @param nbits Integer. Number of bits in each chromosome.
#'
#' @details
#' Bits are sampled independently. Each bit takes the value 0 or 1 with equal
#' probability.
#'
#' @return
#' A numeric matrix with npop rows and nbits columns containing only 0 and 1.
#' Each row corresponds to one chromosome.
#'
#' @author Zhonghui Huang
#'
#' @examples
#' create.pop(npop = 10, nbits = 12)
#'
#' @export
create.pop <- function(npop,
                       nbits) {
  population <- matrix(as.double(NA),
                       nrow = npop,
                       ncol = nbits)

  for (j in 1:npop) {
    population[j,] <- sample(
      x = c(0, 1),
      size = nbits,
      prob = c(0.5, 0.5),
      replace = TRUE
    )
  }
  return(population)
}

#' Tournament selection
#'
#' Select individuals for the next generation using tournament selection.
#'
#' @param data.pop A data.frame containing the current population. The first
#'   nbits columns must be the chromosome (typically coded as 0/1). The data
#'   frame must also contain a column named \code{rank}, where smaller values
#'   indicate better individuals (e.g., rank 1 is best).
#' @param npop Integer. Number of individuals to select for the next generation.
#' @param nbits Integer. Number of bits (genes) per chromosome; i.e., the number
#'   of columns taken from \code{data.pop} to form the chromosome matrix.
#'
#' @return A matrix of dimension npop x nbits containing the selected
#'   chromosomes for the next generation.
#'
#' @author Zhonghui Huang
#'
#' @examples
#' data.pop <- data.frame(fitness = stats::runif(10), rank = rank(stats::runif(10)))
#' population <- matrix(sample(0:1, 100, replace = TRUE), nrow = 10)
#' ga.sel.tournament(data.pop=cbind(as.data.frame(population), data.pop), npop=10, nbits=10)
#'
#' @export
ga.sel.tournament <- function(data.pop, npop, nbits) {
  population <- as.matrix(data.pop[, 1:nbits])
  sel.population <- NULL
  for (tour.cycle in 1:npop) {
    test1 <- data.pop[tour.cycle,]
    candidates <- setdiff(1:nrow(data.pop), tour.cycle)
    test2 <- sample(candidates, size = 1)
    test2.v <- data.pop[test2,]
    if (test2.v$rank <= test1$rank) {
      sel.population.s <- population[test2, , drop = FALSE]
    } else {
      sel.population.s <- population[tour.cycle, , drop = FALSE]
    }
    sel.population <- rbind(sel.population, sel.population.s)
  }
  return(sel.population)
}

#' Crossover operator (one- or two-point) for binary chromosomes
#'
#' Apply one- or two-point crossover to a selected population of binary chromosomes.
#'
#' @param sel.population Numeric matrix of dimension npop by nbits. Each row is a
#'   chromosome and is expected to contain binary values (0/1).
#' @param pcross Single numeric value in \eqn{[0, 1]} giving the probability of applying
#'   crossover to each parent pair.
#' @param npop Single positive even integer giving the population size.
#' @param nbits Single positive integer giving the chromosome length.
#'
#' @details
#' Parents are paired sequentially (1 with 2, 3 with 4, etc.). For each pair,
#' crossover is applied with probability \code{pcross}; otherwise the parents are copied
#' unchanged. Crossover points are sampled from 0:\code{nbits}. If the sampled points do
#' not yield a valid crossover, no crossover is performed for that pair.
#'
#' @return Numeric matrix containing the children population after crossover.
#'
#' @author Zhonghui Huang
#'
#' @examples
#' sel.population <- matrix(sample(0:1, 100, replace = TRUE), nrow = 10)
#' ga.crossover(sel.population = sel.population, pcross = 0.7, npop = 10, nbits = 10)
#'
#' @export

ga.crossover <- function(sel.population,
                         pcross,
                         npop,
                         nbits) {
  children.cross <- NULL
  target <- seq(1, npop, 2)
  for (chromo in target) {
    # Create a blank pair of children chromosome
    children <- matrix(as.double(NA), nrow = 2, ncol = nbits)
    # Select the target pair of parent chromosome.
    parents.cross = sel.population[chromo:(chromo + 1),]
    # Test whether need to do crossover based on the pre-defined probability
    cross.q <- stats::runif(1, min = 0, max = 1)
    if (cross.q > pcross) {
      # No crossover
      children <- parents.cross
    }
    # Crossover
    else{
      # Random select the 2 crossover points for the pair of chromosome.
      cross.points <- sample(0:nbits, size = 2, replace = T)
      cross.point.1 <- min(cross.points)
      cross.point.2 <- max(cross.points)
      # Only one valid crossover point
      if (cross.point.1 == cross.point.2 ||
          cross.point.2 == nbits || cross.point.1 == 0) {
        # no valid crossover point
        if (cross.point.1 == 0 & cross.point.2 == nbits) {
          children[1:2,] = parents.cross
        }
        else if (cross.point.1 == 0 & cross.point.2 == 0) {
          children[1:2,] = parents.cross
        }
        else if (cross.point.1 == nbits & cross.point.2 == nbits) {
          children[1:2,] = parents.cross
        }
        # only one valid crossover point
        else{
          if (cross.point.1 == 0) {
            cross.point = cross.point.2
          }
          else if (cross.point.2 == nbits) {
            cross.point = cross.point.1
          }
          else {
            cross.point = cross.point.1   #e.g., cross.point.1= cross.point.2= 3
          }
          children[1,] <- c(parents.cross[1, 1:cross.point],
                            parents.cross[2, (cross.point + 1):nbits])
          children[2,] <- c(parents.cross[2, 1:cross.point],
                            parents.cross[1, (cross.point + 1):nbits])
        }
      }
      # Two valid crossover points
      else{
        {
          children[1,] <- c(parents.cross[1, 1:cross.point.1],
                            parents.cross[2, (cross.point.1 + 1):cross.point.2],
                            parents.cross[1, (cross.point.2 + 1):nbits])


          children[2,] <- c(parents.cross[2, 1:cross.point.1],
                            parents.cross[1, (cross.point.1 + 1):cross.point.2],
                            parents.cross[2, (cross.point.2 + 1):nbits])
        }
      }
    }
    children.cross <- rbind(children.cross, children)
  }
  return(children.cross)
}


#' Mutation operator for binary genetic algorithms
#'
#' Mutate a binary population by flipping bits with probability pmut.
#'
#' @param children.cross Numeric matrix containing the child population. Rows are
#'   individuals and columns are bits. Values are expected to be 0/1.
#' @param pmut Single numeric value in \eqn{[0, 1]} giving the per-bit mutation probability.
#'
#' @details
#' Mutation is applied independently to each bit (gene). For each position, a
#' Bernoulli trial with success probability pmut determines whether the bit is
#' flipped (0 becomes 1, 1 becomes 0).
#'
#' @return Numeric matrix containing the mutated population.
#'
#' @author Zhonghui Huang
#'
#' @examples
#' children.cross <- matrix(sample(0:1, 120, replace = TRUE), nrow = 10)
#' ga.mutation(children.cross, pmut = 0.1)
#'
#' @export

ga.mutation <- function(children.cross, pmut) {
  # Create a logical mask: TRUE where mutation occurs
  mutation.mask <- matrix(
    stats::rbinom(length(children.cross), 1, pmut) == 1,
    nrow = nrow(children.cross),
    ncol = ncol(children.cross)
  )
  # Flip bits in positions selected for mutation
  children.cross[mutation.mask] <- 1 - children.cross[mutation.mask]
  return(children.cross)
}
