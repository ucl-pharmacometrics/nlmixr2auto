#' Create ants for each travel
#'
#' Create a specified number of ants and initializes their parameters based on the given conditions and probabilities.
#' If `initialize` is `TRUE`, the ants are generated based on a fixed initial design model. Otherwise, they are sampled based on the current probabilities of different paths.
#'
#'
#' @param no.ants Integer. The number of ants to be created.
#' @param initialize Logical. Whether to initialize the ants with predefined values or sample them based on probabilities.
#' @param node.list Data frame. A list of nodes with probabilities used for sampling parameters.
#'
#' @return Data frame. A data frame with initialized parameters for all ants.
#' @examples
#' \dontrun{
#'   initial.phi <- 1
#'   no.nodes <- 20
#'   node.list.all <- data.frame(
#'     travel = 0,
#'     node.no = seq(1, no.nodes, 1),
#'     local.node.no = c(seq(1, 3, 1),
#'                       seq(0, 1, 1),
#'                       seq(0, 1, 1),
#'                       seq(0, 1, 1),
#'                       seq(0, 1, 1),
#'                       seq(0, 1, 1),
#'                       seq(0, 1, 1),
#'                       seq(0, 1, 1),
#'                       seq(1, 3, 1)),
#'     node.names = c("1Cmpt", "2Cmpt", "3Cmpt",
#'                    "eta.vp2.no", "eta.vp2.yes",
#'                    "eta.q2.no", "eta.q2.yes",
#'                    "eta.vp.no", "eta.vp.yes",
#'                    "eta.q.no", "eta.q.yes",
#'                    "eta.vc.no", "eta.vc.yes",
#'                    "mm.no", "mm.yes",
#'                    "mcorr.no", "mcorr.yes",
#'                    "add", "prop", "comb"),
#'     node.group = c(rep(1, 3),
#'                    sort(rep(seq(2, 8, 1), 2)),
#'                    rep(9, 3)),
#'     phi = rep(initial.phi, no.nodes),
#'     delta_phi = rep(0, no.nodes),
#'     p = c(rep(round(1 / 3, 3), 3),
#'           rep(0.5, 14),
#'           rep(round(1 / 3, 3), 3))
#'   )
#'   ants <- create.ant(no.ants = 10, initialize = TRUE, node.list = node.list.all)
#' }
#' @export


create.ant <- function(search.space,
                       no.ants,
                       initialize,
                       node.list) {
  if (search.space == 1) {
    no.nodes = 22
    no.ants <- no.ants
    
    cmpt.iv.r1      <- data.frame(matrix(nrow = 1, ncol = no.ants))
    eta.km.r1       <- data.frame(matrix(nrow = 1, ncol = no.ants))
    eta.vc.r1       <- data.frame(matrix(nrow = 1, ncol = no.ants))
    eta.vp.r1       <- data.frame(matrix(nrow = 1, ncol = no.ants))
    eta.vp2.r1      <- data.frame(matrix(nrow = 1, ncol = no.ants))
    eta.q.r1        <- data.frame(matrix(nrow = 1, ncol = no.ants))
    eta.q2.r1       <- data.frame(matrix(nrow = 1, ncol = no.ants))
    mm.r1           <- data.frame(matrix(nrow = 1, ncol = no.ants))
    eta.mcorr.r1    <- data.frame(matrix(nrow = 1, ncol = no.ants))
    eta.rv.r1       <- data.frame(matrix(nrow = 1, ncol = no.ants))
    
    
    if (initialize == T) {
      # Start from simple model
      cmpt.iv.r1[1:6] <- c(1, 1, 2, 2, 3, 3)
      eta.vc.r1[1:6] <- c(0, 0, 0, 0, 0, 0)
      eta.vp.r1[1:6] <- c(-1,-1, 0, 0, 0, 0)
      eta.q.r1[1:6] <- c(-1,-1, 0, 0, 0, 0)
      eta.vp2.r1[1:6] <- c(-1,-1,-1,-1, 0, 0)
      eta.q2.r1[1:6] <- c(-1,-1,-1,-1, 0, 0)
      eta.rv.r1[1:6] <- c(1, 2, 1, 2, 1, 2)
      eta.mcorr.r1[1:6] <- c(0, 0, 0, 0, 0, 0)
      mm.r1[1:6] <- c(0, 0, 0, 0, 0, 0)
      eta.km.r1[1:6] <- c(-1, -1, -1, -1, -1, -1)
      
      if (no.ants > 6) {
        cmpt.iv.r1[7:no.ants] <- sample(
          x = seq(1, 3, 1),
          size = (no.ants - 6),
          prob = c(node.list[1:3,]$p),
          replace = T
        )
        
        eta.vc.r1[7:no.ants] <- sample(
          x = c(0, 1),
          size = (no.ants - 6),
          prob = c(node.list[12:13,]$p),
          replace = T
        )
        
        mm.r1[7:no.ants] <- sample(
          x = c(0, 1),
          size = (no.ants - 6),
          prob = c(node.list[14:15,]$p),
          replace = T
        )
        
        eta.mcorr.r1[7:no.ants] <- sample(
          x = c(0, 1),
          size = (no.ants - 6),
          prob = c(node.list[18:19,]$p),
          replace = T
        )
        
        eta.rv.r1[7:no.ants] <- sample(
          x = seq(1, 3, 1),
          size = (no.ants - 6),
          prob = c(node.list[20:22,]$p),
          replace = T
        )
        
        # identify n/a situation, use -1
        
        for (j in 7:no.ants) {
          if (cmpt.iv.r1[j] == 3) {
            eta.vp2.r1[j] <- sample(
              x = c(0, 1),
              size = 1,
              prob = c(node.list[4:5,]$p)
            )
            eta.q2.r1[j] <-
              sample(
                x = c(0, 1),
                size = 1,
                prob = c(node.list[6:7,]$p)
              )
            
          }
          else{
            eta.vp2.r1[j] <- -1
            eta.q2.r1[j] <- -1
          }
        }
        
        
        for (j in 7:no.ants) {
          if (cmpt.iv.r1[j] > 1) {
            eta.vp.r1[j] <- sample(
              x = c(0, 1),
              size = 1,
              prob = c(node.list[8:9,]$p)
            )
            eta.q.r1[j] <-
              sample(
                x = c(0, 1),
                size = 1,
                prob = c(node.list[10:11,]$p)
              )
          }
          else{
            eta.vp.r1[j] <- -1
            eta.q.r1[j] <- -1
          }
        }
        
        
        for (j in 7:no.ants) {
          
          if (mm.r1[j] ==1) {
           eta.km.r1[j] <- sample(
           x = c(0, 1),
           size = 1,
           prob = c(node.list[16:17,]$p),
           replace = T
           )
          }
          else{
            eta.km.r1[j] <- -1
            
          }
    
       }
      
      } # close else
      
    } # close initialize
    
    else{
      cmpt.iv.r1[1:no.ants] <- sample(
        x = seq(1, 3, 1),
        size = no.ants,
        prob = c(node.list[1:3,]$p),
        replace = T
      )
      
      eta.vc.r1[1:no.ants] <- sample(
        x = c(0, 1),
        size = no.ants,
        prob = c(node.list[12:13,]$p),
        replace = T
      )
      
      mm.r1[1:no.ants] <- sample(
        x = c(0, 1),
        size = no.ants,
        prob = c(node.list[14:15,]$p),
        replace = T
      )
      
      eta.km.r1[1:no.ants] <- sample(
        x = c(0, 1),
        size = no.ants,
        prob = c(node.list[16:17,]$p),
        replace = T
      )
      
      eta.mcorr.r1[1:no.ants] <- sample(
        x = c(0, 1),
        size = no.ants,
        prob = c(node.list[18:19,]$p),
        replace = T
      )
      
      eta.rv.r1[1:no.ants] <- sample(
        x = seq(1, 3, 1),
        size = no.ants,
        prob = c(node.list[20:22,]$p),
        replace = T
      )
      
      
      
      for (j in 1:no.ants) {
        if (cmpt.iv.r1[j] == 3) {
          eta.vp2.r1[j] <- sample(
            x = c(0, 1),
            size = 1,
            prob = c(node.list[4:5,]$p)
          )
          eta.q2.r1[j] <-
            sample(
              x = c(0, 1),
              size = 1,
              prob = c(node.list[6:7,]$p)
            )
          
        }
        else{
          eta.vp2.r1[j] <- -1
          eta.q2.r1[j] <- -1
        }
      }
      
      
      # 2cmpt, 3cmpt
      for (j in 1:no.ants) {
        if (cmpt.iv.r1[j] > 1) {
          eta.vp.r1[j] <- sample(
            x = c(0, 1),
            size = 1,
            prob = c(node.list[8:9,]$p)
          )
          eta.q.r1[j] <-
            sample(
              x = c(0, 1),
              size = 1,
              prob = c(node.list[10:11,]$p)
            )
        }
        else{
          eta.vp.r1[j] <- -1
          eta.q.r1[j] <- -1
        }
        
      }
      
      
      # 2cmpt, 3cmpt
      for (j in 1:no.ants) {
        if (mm.r1[j] == 1) {
          eta.km.r1[j] <- sample(
            x = c(0, 1),
            size = 1,
            prob = c(node.list[16:17,]$p)
          )
    
        }
        else{
          eta.km.r1[j] <- -1
      
        }
      }
          
    }
    
  }# close if search.space<-1
  

  ants.all <- rbind(
    cmpt.iv.r1,
    eta.vp2.r1,
    eta.q2.r1,
    eta.vp.r1,
    eta.q.r1,
    eta.vc.r1,
    mm.r1,
    eta.km.r1,
    eta.mcorr.r1,
    eta.rv.r1
  )
  
  return(ants.all)
  
}
