% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runlocal.R
\name{runlocal}
\alias{runlocal}
\title{Perform local exhaustive search for genetic algorithm}
\usage{
runlocal(
  ga.iter,
  dat,
  search.space,
  sel.best.code,
  nbits,
  sig.diff,
  crse,
  cshrink,
  lbcl,
  lbvc,
  lbvp,
  lbq,
  lbvp2,
  lbq2,
  cadd,
  cprop,
  ccorr,
  penalty.type,
  penalty.value,
  ...
)
}
\arguments{
\item{ga.iter}{Current genetic algorithm iteration number.}

\item{dat}{Data frame containing the input data for model fitting.}

\item{search.space}{An integer indicating the search space.}

\item{sel.best.code}{The best chromosome from the previous iteration.}

\item{nbits}{Number of bits (parameters) in the chromosome.}

\item{sig.diff}{Significant difference threshold for ranking.}

\item{crse}{Constraint value for relative standard error.}

\item{cshrink}{Constraint value for shrinkage.}

\item{lbcl}{Lower bound for clearance (CL).}

\item{lbvc}{Lower bound for volume of central compartment (Vc).}

\item{lbvp}{Lower bound for volume of peripheral compartment (Vp).}

\item{lbq}{Lower bound for inter-compartmental clearance (Q).}

\item{lbvp2}{Lower bound for volume of second peripheral compartment (Vp2).}

\item{lbq2}{Lower bound for second inter-compartmental clearance (Q2).}

\item{cadd}{Constraint value for additive error.}

\item{cprop}{Constraint value for proportional error.}

\item{ccorr}{Constraint value for correlation.}

\item{penalty.type}{Type of penalty applied during fitness evaluation.}

\item{penalty.value}{Value of the penalty applied during fitness evaluation.}

\item{...}{Additional arguments passed to \code{ga.mod.run}.}

\item{modi}{Current model iteration number.}
}
\value{
A data frame containing the population of chromosomes generated during the local exhaustive search and their fitness values.
}
\description{
Generates a population of chromosomes for local exhaustive search based on the best individual
from the previous iteration, and valuates the fitness of each new chromosome and ranks them.
}
\examples{
\dontrun{

dat<-pheno_sd
ga.iter<-1
search.space<-1
sel.best.code <- data.frame(0,
                           0,
                            0,
                           0,
                           0,
                           0,
                           0,
                           0,
                           0,
                           0,
                           1,
                           1)
sig.diff<-1
nbits<-12
crse <- 20
cshrink <- 20
lbcl <- 0.1
lbvc <- 0.1
lbvp <- 0.1
lbq <- 0.1
lbvp2 <- 0.1
lbq2 <- 0.1
cadd <- 0.1
cprop <- 0.1
ccorr <- 0.1
penalty.type <- 3
penalty.value <- 10000

 ls.population<-runlocal( ga.iter,
                       dat,
                       search.space,
                       sel.best.code,
                        nbits,
                       sig.diff,
                       crse,
                       cshrink,
                       lbcl,
                       lbvc,
                       lbvp,
                       lbq,
                       lbvp2,
                       lbq2,
                       cadd,
                       cprop,
                        ccorr,
                       penalty.type,
                       penalty.value,
                        control = saemControl(
                         seed = 1234,
                         print = 5,
                         nBurn = 20,
                         nEm = 30,
                         logLik = T,
                         rxControl = rxControl(cores = 2,
                                               maxsteps =70000)))

}

}
