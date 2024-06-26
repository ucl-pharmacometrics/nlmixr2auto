#' Tournament selection for genetic algorithm
#'
#' Perform tournament selection to choose individuals for the next generation
#' in a genetic algorithm.
#'
#' @param data.pop A data frame containing the current population with their fitness values and rankings.
#' @param population A matrix representing the current population of chromosomes.
#' @param npopsize An integer specifying the size of the population.
#' @param nbits An integer specifying the number of bits in each chromosome.
#'
#' @return A matrix representing the selected population for the next generation.
#'
#' @examples

#' data.pop <- data.frame(fitness = runif(10), rank = rank(runif(10)))
#' population <- matrix(sample(0:1, 100, replace = TRUE), nrow = 10)
#' npopsize <- 10
#' nbits <- 10
#'ga.sel.tournament(data.pop, population, npopsize, nbits)

#'
#' @export
ga.sel.tournament <- function(data.pop,
                              population,
                              npopsize,
                              nbits) {
  sel.population <- NULL
  for (tour.cycle in 1:(npopsize)) {
    # Select first indvidual
    if (tour.cycle <= nrow(data.pop)) {
      test1 = data.pop[tour.cycle,]
    }
    # Randomly select the second individual, but should be different with the first one
    for (ntest2 in 1:10000) {
      test2 = sample(1:nrow(data.pop), size = 1)
      if (test2 == tour.cycle) {
        test2 = sample(1:nrow(data.pop), size = 1)
      }
      else{
        break
      }
    }
    # Get the fitness information
    test2.v <- data.pop[test2,]
    # Compare the rankings between two individuals, one with lower rankings is to survive.
    if (test2.v$rank <= test1$rank) {
      sel.population.s <- population[test2, , drop = F]
    }
    if (test2.v$rank > test1$rank) {
      sel.population.s <-
        population[as.numeric(rownames(test1)), , drop = F]
    }
    
    sel.population <- rbind(sel.population, sel.population.s)
  }
  return(sel.population)
}

#' Crossover in genetic algorithm
#'
#' Perform a crossover operation on a selected population in a genetic algorithm.
#' Use a predefined probability to decide whether to perform crossover on pairs of parent chromosomes
#' and generates children chromosomes by combining segments of the parent chromosomes.
#'
#' @param sel.population A matrix representing the selected population of chromosomes.
#' @param prob.crossover The probability of crossover in the genetic algorithm.
#' @param npopsize An integer specifying the size of the population.
#' @param nbits An integer specifying the number of bits in each chromosome.
#' @param bit.names A vector of character strings representing the names of the bits.
#'
#' @return A matrix representing the population of children chromosomes after crossover.
#'
#' @examples

#' sel.population <- matrix(sample(0:1, 100, replace = TRUE), nrow = 10)
#' prob.crossover <- 0.8
#' npopsize <- 10
#' nbits <- 12
#  bit.names <- c("cmpt.iv1",
#                "cmpt.iv2",
#                "eta.km",
#                "eta.vc",
#                "eta.vp",
#                "eta.vp2",
#                "eta.q",
#                "eta.q2",
#                "mm",
#                "mcorr",
#                "rv1",
#                "rv2")
#' ga.crossover(sel.population, prob.crossover, npopsize, nbits, bit.names)
#'
#' @export

ga.crossover <- function(sel.population,
                         prob.crossover,
                         npopsize,
                         nbits,
                         bit.names) {
  children.all <- NULL
  # Sequential pairing
  target <- seq(1, npopsize, 2)
  for (chromo in target) {
    # Create a blank pair of children chromosome
    children <- matrix(as.double(NA), nrow = 2, ncol = nbits)
    # Select the target pair of parent chromosome.
    parents.cross = sel.population[chromo:(chromo + 1),]
    # Test whether need to do crossover based on the pre-defined probability
    cross.q <- runif(1, min = 0, max = 1)
    if (cross.q > prob.crossover) {
      # No crossover
      children <- parents.cross
    }
    # Conduct crossover
    else{
      # Random select the 2 crossover points for the pair of chromosome.
      # This method includes the scenario of only one crossover point selected
      cross.points <- sample(0:nbits, size = 2, replace = T)
      # print(cross.points)
      # print(parents.cross)
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
    children.all <- rbind(children.all, children)
  }
  colnames(children.all) <- bit.names
  return(children.all)
}

#' Mutation in Genetic Algorithm
#'
#' Perform a mutation operation on a population of chromosomes in a genetic algorithm.
#' Each bit in the chromosome has a chance to mutate based on a predefined mutation probability.
#'
#' @param children.all A matrix representing the population of children chromosomes.
#' @param prob.mutation The probability of mutation for each bit in the chromosome.
#' @param npopsize An integer specifying the size of the population.
#' @param nbits An integer specifying the number of bits in each chromosome.
#'
#' @return A matrix representing the population of children chromosomes after mutation.
#'
#' @examples
#' \dontrun{
#' children.all <- matrix(sample(0:1, 120, replace = TRUE), nrow = 10)
#' prob.mutation <- 0.1
#' npopsize <- 10
#' nbits <- 12
#' mutated_population <- ga.mutation(children.all, prob.mutation, npopsize, nbits)
#' print(mutated_population)
#' }
#'
#' @export

ga.mutation <- function(children.all,
                        prob.mutation,
                        npopsize,
                        nbits) {
  for (mutate.j in 1:nbits) {
    for (mutate.i in 1:npopsize) {
      # 0 means no mutation, 1 means mutation
      prob.test = sample(
        x = c(0, 1),
        size = 1,
        replace = T,
        prob = c((1 - prob.mutation), prob.mutation)
      )
      if (prob.test == 0) {
        children.all[mutate.i , mutate.j] = children.all[mutate.i, mutate.j]
      }
      else{
        if (children.all[mutate.i , mutate.j] == 0) {
          children.all[mutate.i , mutate.j] = 1
        }
        else{
          children.all[mutate.i, mutate.j] = 0
        }
      }
    }
  }
  return(children.all)
}
