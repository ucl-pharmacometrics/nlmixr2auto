#' Calculate pheromone update for ants
#'
#' Calculate the updated pheromone values for each node based on the results of the current iteration and the pheromone evaporation rate.
#'
#' @param r An integer representing the current iteration or travel.
#' @param search.space An integer specifying the search space type (1 for standard).
#' @param input.ants.dat A data frame containing the data of the input ants.
#' @param all.lib.dat A data frame containing the library data with fitness values.
#' @param node.list.0 A template data frame containing the all node pheromone information.
#' @param node.list.all A data frame containing the all node pheromone information.
#' @param alpha.value A numeric value for the alpha parameter in the pheromone calculation.
#' @param rho A numeric value for the evaporation rate of pheromone.
#' @param sig.diff A numeric value used in the ranking calculation.
#' @param lower.limit.phi A numeric value specifying the lower limit for the pheromone level.
#' @param upper.limit.phi A numeric value specifying the upper limit for the pheromone level.
#'
#' @return Data frame. The updated node list with new pheromone values.
#' @examples
#' \dontrun{
#' r <- 1
#' search.space <- 1
#' no.ants <- 10
#' node.list.0 <- data.frame()  # Define node.list.0 as per your application needs
#' input.ants.dat <- create.ant(search.space = search.space,
#'                              no.ants = no.ants,
#'                              initialize = TRUE,
#'                              node.list = node.list.0)
#'
#' all.lib.dat <- data.frame(fitness = runif(10)) # Example library data
#' no.nodes <- 22
#' initial.phi <- 1
#' node.list.all <- data.frame(
#'   travel = 0,
#'   node.no = seq(1, no.nodes, 1),
#'   local.node.no = c(seq(1, 3, 1), seq(0, 1, 1), seq(0, 1, 1), seq(0, 1, 1), seq(0, 1, 1),
#'                     seq(0, 1, 1), seq(0, 1, 1), seq(0, 1, 1), seq(0, 1, 1), seq(1, 3, 1)),
#'   node.names = c("1Cmpt", "2Cmpt", "3Cmpt", "eta.vp2.no", "eta.vp2.yes", "eta.q2.no", "eta.q2.yes",
#'                  "eta.vp.no", "eta.vp.yes", "eta.q.no", "eta.q.yes", "eta.vc.no", "eta.vc.yes",
#'                  "mm.no", "mm.yes", "eta.km.no", "eta.km.yes", "mcorr.no", "mcorr.yes",
#'                  "add", "prop", "comb"),
#'   node.group = c(rep(1, 3), sort(rep(seq(2, 9, 1), 2)), rep(10, 3)),
#'   phi = rep(initial.phi, no.nodes),
#'   delta_phi = rep(0, no.nodes),
#'   p = c(rep(round(1 / 3, 3), 3), rep(0.5, 16), rep(round(1 / 3, 3), 3))
#' )
#' node.list.0=node.list.all
#' alpha.value <- 1
#' rho <- 0.2
#' sig.diff <- 0.1
#' lower.limit.phi <- 1
#' upper.limit.phi <- 10
#'
#' node.list <- phi.calculate(r, search.space, input.ants.dat, all.lib.dat,
#' node.list.all, alpha.value, rho, sig.diff, lower.limit.phi, upper.limit.phi)
#' }
#'
phi.calculate <- function(r,
                          search.space,
                          input.ants.dat,
                          all.lib.dat,
                          node.list.0,
                          node.list.all,
                          alpha.value,
                          rho,
                          sig.diff,
                          lower.limit.phi,
                          upper.limit.phi) {
  if (search.space == 1) {
    length.string <- 10
    
    # create the new node.list for finished travel
    node.list <- node.list.0
    node.list$travel = r
    phi.dat <<- input.ants.dat
    
    all.lib.dat$allrank <- rank_new(all.lib.dat$fitness,
                                    sig.diff = sig.diff)
    
    write.csv(all.lib.dat, file = paste0("all.lib.dat.", r, ".csv"))
    
    # extract the pheromone information for the current travel
    row.start <- nrow(all.lib.dat) - ncol(input.ants.dat) + 1
    row.stop <- nrow(all.lib.dat)
    
    
    phi <-
      (1 / (all.lib.dat[row.start:row.stop,]$allrank)) ^ alpha.value
    
    # Replace the number less than 0
    phi.dat[(length.string + 1),] <-
      pmax(phi, 0) # based on the equation, it should be be less than 0
    
    # change the number for single node to the overall node number
    phi.dat[2,] <- phi.dat[2,] + 4
    phi.dat[3,] <- phi.dat[3,] + 6
    phi.dat[4,] <- phi.dat[4,] + 8
    phi.dat[5,] <- phi.dat[5,] + 10
    phi.dat[6,] <- phi.dat[6,] + 12
    phi.dat[7,] <- phi.dat[7,] + 14
    phi.dat[8,] <- phi.dat[8,] + 16
    phi.dat[9,] <- phi.dat[9,] + 18
    phi.dat[10,] <- phi.dat[10,] + 19
    
    # write.csv(phi.dat,file = paste0("phi.dat",r,".csv") )
    
    # Calculate the total remaining pheromone from the previous iteration, then evaporate with the ratio of rho
    for (n in 1:nrow(node.list)) {
      if (n < 4) {
        test.datM <- rbind(phi.dat[1,], phi.dat[11,])
        test.datM[3,] <- test.datM[1,] == n
        node.list[n,]$delta_phi = sum(test.datM[2,] * test.datM[3,])
        # pheromone generation and evaporation
        node.list[n,]$phi = (1 - rho) * sum(node.list.all[node.list.all$node.no ==
                                                            n &
                                                            node.list.all$travel > (r - 2),]$phi) + node.list[n,]$delta_phi
      }
      
      else if (n < 6) {
        test.datM <- rbind(phi.dat[2,], phi.dat[11,])
        test.datM[3,] <- test.datM[1,] == n
        node.list[n,]$delta_phi = sum(test.datM[2,] * test.datM[3,])
        node.list[n,]$phi = (1 - rho) * sum(node.list.all[node.list.all$node.no ==
                                                            n &
                                                            node.list.all$travel > (r - 2),]$phi) + node.list[n,]$delta_phi
      }
      
      else if (n < 8) {
        test.datM <- rbind(phi.dat[3,], phi.dat[11,])
        test.datM[3,] <- test.datM[1,] == n
        node.list[n,]$delta_phi = sum(test.datM[2,] * test.datM[3,])
        node.list[n,]$phi = (1 - rho) * sum(node.list.all[node.list.all$node.no ==
                                                            n &
                                                            node.list.all$travel > (r - 2),]$phi) + node.list[n,]$delta_phi
        
      }
      else if (n < 10) {
        test.datM <- rbind(phi.dat[4,], phi.dat[11,])
        test.datM[3,] <- test.datM[1,] == n
        node.list[n,]$delta_phi = sum(test.datM[2,] * test.datM[3,])
        node.list[n,]$phi = (1 - rho) * sum(node.list.all[node.list.all$node.no ==
                                                            n &
                                                            node.list.all$travel > (r - 2),]$phi) + node.list[n,]$delta_phi
        
      }
      else if (n < 12) {
        test.datM <- rbind(phi.dat[5,], phi.dat[11,])
        test.datM[3,] <- test.datM[1,] == n
        node.list[n,]$delta_phi = sum(test.datM[2,] * test.datM[3,])
        node.list[n,]$phi = (1 - rho) * sum(node.list.all[node.list.all$node.no ==
                                                            n &
                                                            node.list.all$travel > (r - 2),]$phi) + node.list[n,]$delta_phi
        
      }
      else if (n < 14) {
        test.datM <- rbind(phi.dat[6,], phi.dat[11,])
        test.datM[3,] <- test.datM[1,] == n
        node.list[n,]$delta_phi = sum(test.datM[2,] * test.datM[3,])
        node.list[n,]$phi = (1 - rho) * sum(node.list.all[node.list.all$node.no ==
                                                            n &
                                                            node.list.all$travel > (r - 2),]$phi) + node.list[n,]$delta_phi
        
      }
      
      else if (n < 16) {
        test.datM <- rbind(phi.dat[7,], phi.dat[11,])
        test.datM[3,] <- test.datM[1,] == n
        node.list[n,]$delta_phi = sum(test.datM[2,] * test.datM[3,])
        node.list[n,]$phi = (1 - rho) * sum(node.list.all[node.list.all$node.no ==
                                                            n &
                                                            node.list.all$travel > (r - 2),]$phi) + node.list[n,]$delta_phi
      }
      
      else if (n < 18) {
        test.datM <- rbind(phi.dat[8,], phi.dat[11,])
        test.datM[3,] <- test.datM[1,] == n
        node.list[n,]$delta_phi = sum(test.datM[2,] * test.datM[3,])
        node.list[n,]$phi = (1 - rho) * sum(node.list.all[node.list.all$node.no ==
                                                            n &
                                                            node.list.all$travel > (r - 2),]$phi) + node.list[n,]$delta_phi
      }
      
      else if (n < 20) {
        test.datM <- rbind(phi.dat[9,], phi.dat[11,])
        test.datM[3,] <- test.datM[1,] == n
        node.list[n,]$delta_phi = sum(test.datM[2,] * test.datM[3,])
        node.list[n,]$phi = (1 - rho) * sum(node.list.all[node.list.all$node.no ==
                                                            n &
                                                            node.list.all$travel > (r - 2),]$phi) + node.list[n,]$delta_phi
      }
      
      else{
        test.datM <- rbind(phi.dat[10,], phi.dat[11,])
        test.datM[3,] <- test.datM[1,] == n
        node.list[n,]$delta_phi = sum(test.datM[2,] * test.datM[3,])
        node.list[n,]$phi = (1 - rho) * sum(node.list.all[node.list.all$node.no ==
                                                            n &
                                                            node.list.all$travel > (r - 2),]$phi) + node.list[n,]$delta_phi
        
      }
      
      node.list[n,]$phi <- pmax(node.list[n,]$phi, lower.limit.phi)
      node.list[n,]$phi <- pmin(node.list[n,]$phi, upper.limit.phi)
    }
    
  }
  return(node.list)
}




#' Calculate selection probabilities for nodes
#'
#' Calculate the selection probabilities for nodes in the ant colony optimization
#'
#' @param search.space An integer specifying the search space type (1 for standard).
#' @param node.list.cal A data frame containing the list of nodes with their pheromone levels and other attributes.
#' @return A data frame representing the updated node list with calculated selection probabilities.
#' @examples
#' \dontrun{
#' # Example usage:
#' search.space <- 1
#' no.nodes <- 22
#' initial.phi <- 1
#' node.list.cal <- data.frame(
#'   travel = 1,
#'   node.no = seq(1, no.nodes, 1),
#'   local.node.no = c(seq(1, 3, 1), seq(0, 1, 1), seq(0, 1, 1), seq(0, 1, 1), seq(0, 1, 1),
#'                     seq(0, 1, 1), seq(0, 1, 1), seq(0, 1, 1), seq(0, 1, 1), seq(1, 3, 1)),
#'   node.names = c("1Cmpt", "2Cmpt", "3Cmpt", "eta.vp2.no", "eta.vp2.yes", "eta.q2.no", "eta.q2.yes",
#'                  "eta.vp.no", "eta.vp.yes", "eta.q.no", "eta.q.yes", "eta.vc.no", "eta.vc.yes",
#'                  "mm.no", "mm.yes", "eta.km.no", "eta.km.yes", "mcorr.no", "mcorr.yes",
#'                  "add", "prop", "comb"),
#'   node.group = c(rep(1, 3), sort(rep(seq(2, 9, 1), 2)), rep(10, 3)),
#'   phi = seq(1,22),
#'   delta_phi = rep(0, no.nodes),
#'   p = c(rep(round(1 / 3, 3), 3), rep(0.5, 16), rep(round(1 / 3, 3), 3))
#' )
#' updated.node.list.cal <- p.calculation(search.space, node.list.cal)
#' }
#'

p.calculation <- function(search.space,
                          node.list.cal) {
  if (search.space == 1) {
    for (n in 1:nrow(node.list.cal)) {
      node.list.cal[n,]$p <-
        node.list.cal[n,]$phi / sum(node.list.cal[node.list.cal$node.group == node.list.cal[n,]$node.group,]$phi)
      # Ensure a 20% chance of selection.
      # for the node group with two local node options
      if (node.list.cal[n,]$node.group %in% c(2:9)) {
        node.list.cal[n,]$p <- pmax(node.list.cal[n,]$p, 0.2)
        node.list.cal[n,]$p <- pmin(node.list.cal[n,]$p, 0.8)
      }
      
      # for the node group with three local node options
      if (node.list.cal[n,]$node.group %in% c(1, 10)) {
        node.list.cal[n,]$p <- pmax(node.list.cal[n,]$p, 0.25)
        node.list.cal[n,]$p <- pmin(node.list.cal[n,]$p, 0.75)
        
      }
      
    }
    
    return(node.list.cal)
    
  }
}
