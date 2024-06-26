% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create.ant.R
\name{create.ant}
\alias{create.ant}
\title{Create ants for each travel}
\usage{
create.ant(search.space, no.ants, initialize, node.list)
}
\arguments{
\item{no.ants}{Integer. The number of ants to be created.}

\item{initialize}{Logical. Whether to initialize the ants with predefined values or sample them based on probabilities.}

\item{node.list}{Data frame. A list of nodes with probabilities used for sampling parameters.}
}
\value{
Data frame. A data frame with initialized parameters for all ants.
}
\description{
Create a specified number of ants and initializes their parameters based on the given conditions and probabilities.
If \code{initialize} is \code{TRUE}, the ants are generated based on a fixed initial design model. Otherwise, they are sampled based on the current probabilities of different paths.
}
\examples{
\dontrun{
  initial.phi <- 1
  no.nodes <- 20
  node.list.all <- data.frame(
    travel = 0,
    node.no = seq(1, no.nodes, 1),
    local.node.no = c(seq(1, 3, 1),
                      seq(0, 1, 1),
                      seq(0, 1, 1),
                      seq(0, 1, 1),
                      seq(0, 1, 1),
                      seq(0, 1, 1),
                      seq(0, 1, 1),
                      seq(0, 1, 1),
                      seq(1, 3, 1)),
    node.names = c("1Cmpt", "2Cmpt", "3Cmpt",
                   "eta.vp2.no", "eta.vp2.yes",
                   "eta.q2.no", "eta.q2.yes",
                   "eta.vp.no", "eta.vp.yes",
                   "eta.q.no", "eta.q.yes",
                   "eta.vc.no", "eta.vc.yes",
                   "mm.no", "mm.yes",
                   "mcorr.no", "mcorr.yes",
                   "add", "prop", "comb"),
    node.group = c(rep(1, 3),
                   sort(rep(seq(2, 8, 1), 2)),
                   rep(9, 3)),
    phi = rep(initial.phi, no.nodes),
    delta_phi = rep(0, no.nodes),
    p = c(rep(round(1 / 3, 3), 3),
          rep(0.5, 14),
          rep(round(1 / 3, 3), 3))
  )
  ants <- create.ant(no.ants = 10, initialize = TRUE, node.list = node.list.all)
}
}
