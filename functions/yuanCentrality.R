##
## wdnet: Weighted directed network
## Copyright (C) 2021  Panpan Zhang and Jun Yan
## Jun Yan <jun.yan@uconn.edu>
##
## This file is part of the R package wdnet.
##
## The R package wdnet is free software: You can redistribute it and/or
## modify it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or any later
## version (at your option). See the GNU General Public License at
## <https://www.gnu.org/licenses/> for details.
##
## The R package wdnet is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##

#' @importFrom igraph distances graph_from_adjacency_matrix
#' @importFrom rARPACK eigs
NULL

#' Degree-based centrality
#'
#' Compute the degree centrality measures of the vertices in a weighted and directed
#' network represented through its adjacency matrix.
#'
#' @param adj is an adjacency matrix of an weighted and directed network
#' @param alpha is a tuning parameter. The value of alpha must be nonnegative. By convetion, 
#' alpha takes a value from 0 to 1 (default).
#' @param mode which mode to compute: "out" (default) or "in"? For undirected networks, this
#' setting is irrelavent.
#'
#' @return a list of node names and associated degree centrality measures
#'
#' @references
#' \itemize{
#' \item Opsahl, T., Agneessens, F., Skvoretz, J. (2010). Node centrality 
#' in weighted networks: Generalizing degree and shortest paths. 
#' \emph{Social Networks}, 32, 245--251.
#' \item Zhang, P., Zhao, J. and Yan, J. (2020+) Centrality measures of 
#' networks with application to world input-output tables
#' }
#'
#' @note 
#' Function \code{degree_c} is an extension of function \code{strength} 
#' in package \code{igraph} and an altertanive of function \code{degree_w} in 
#' package \code{tnet}. Function \code{degree_c} uses adjacency matrix as
#' input.
#'
#' @examples
#' ## Generate a network according to the Erd\"{o}s-Renyi model of order 20
#' ## and parameter p = 0.3
#' edge_ER <- rbinom(400,1,0.3)
#' weight_ER <- sapply(edge_ER, function(x) x*sample(3,1))
#' adj_ER <- matrix(weight_ER,20,20)
#' mydegree <- degree_c(adj_ER, alpha = 0.8, mode = "in")
#' system.time(mydegree)
#' 
#' @export

degree_c <- function(adj, alpha = 1, mode = "out"){
  if (alpha < 0){
    stop("The tuning parameter alpha must be nonnegative!")
  }
  if (dim(adj)[1]!=dim(adj)[2]){
    stop("The adjacency matrix must be a square matrix!")
  }
  else{
    if (isSymmetric(adj) == TRUE){
      warnings("The analyzed network is undirected!")
    }
    deg_c_output <- matrix(NA, nrow = dim(adj)[1], ncol = 2)
    adj_name <- colnames(adj)
    if (is.null(adj_name) == FALSE){
      deg_c_output <- adj_name
    }
    else{
      deg_c_output[,1] <- c(1:dim(adj)[1])
    }
    colnames(deg_c_output) <- c("name","deg_c")
    adj_deg <- adj
    adj_deg[which(adj_deg > 0)] <- 1
    if (mode == "in"){
      deg_c_output[,2] <- colSums(adj)^alpha + colSums(adj_deg)^(1 - alpha)
    }
    if (mode == "out"){
      deg_c_output[,2] <- rowSums(adj)^alpha + rowSums(adj_deg)^(1 - alpha)
    }
    return(deg_c_output)
  }
}

#' Closeness centrality
#'
#' Compute the closeness centrality measures of the vertices in a weighted and directed 
#' network represented through its adjacency matrix.
#' 
#' @param adj is an adjacency matrix of an weighted and directed network
#' @param alpha is a tuning parameter. The value of alpha must be nonnegative. By convetion, 
#' alpha takes a value from 0 to 1 (default).
#' @param mode which mode to compute: "out" (default) or "in"? For undirected networks, this
#' setting is irrelavent.
#' @param method which method to use: "harmonic" (default) or "standard"?
#' @param distance whether to consider the entries in the adjacency matrix as distances or
#' strong connections. The default setting is \code{FALSE}.
#' 
#' @return a list of node names and associated closeness centrality measures
#'
#' @references
#' \itemize{
#' \item Dijkstra, E.W. (1959). A note on two problems in connexion with 
#' graphs. \emph{Numerische Mathematik}, 1, 269--271.
#' \item Newman, M.E.J. (2003). The structure and function of complex
#' networks. \emph{SIAM review}, 45(2), 167--256.
#' \item Opsahl, T., Agneessens, F., Skvoretz, J. (2010). Node centrality 
#' in weighted networks: Generalizing degree and shortest paths. 
#' \emph{Social Networks}, 32, 245--251.
#' \item Zhang, P., Zhao, J. and Yan, J. (2020+) Centrality measures of 
#' networks with application to world input-output tables
#' }
#'
#' @note 
#' Function \code{closeness_c} is an extension of function \code{closeness} 
#' in package \code{igraph} and function \code{closeness_w} in 
#' package \code{tnet}. The method of computing distances between verticies
#' is the \emph{Dijkstra's algorithm}.
#'
#' @examples
#' ## Generate a network according to the Erd\"{o}s-Renyi model of order 20
#' ## and parameter p = 0.3
#' edge_ER <- rbinom(400,1,0.3)
#' weight_ER <- sapply(edge_ER, function(x){x*sample(3,1)})
#' adj_ER <- matrix(weight_ER,20,20)
#' myclose <- closeness_c(adj_ER, alpha = 0.8, mode = "out",
#' method = "harmonic", distance = FALSE)
#' system.time(myclose)
#' 
#' @export

closeness_c <- function(adj, alpha = 1, mode = "out",
                        method = "harmonic", distance = FALSE){
  if (alpha < 0){
    stop("The tuning parameter alpha must be nonnegative!")
  }
  if (dim(adj)[1]!=dim(adj)[2]){
    stop("The adjacency matrix must be a square matrix!")
  }
  else{
    closeness_c_output <- matrix(NA, nrow = dim(adj)[1], ncol = 2)
    adj_name <- colnames(adj)
    if (is.null(adj_name) == FALSE){
      closeness_c_output[,1] <- adj_name
    }
    else{
      closeness_c_output[,1] <- c(1:dim(adj)[1])
    }
    colnames(closeness_c_output) <- c("name","closeness_c")
    if (distance == FALSE){
      adj <- (1/adj)^alpha
    } else if (distance == TRUE){
      adj <- adj^alpha
    }
    temp_g <- igraph::graph_from_adjacency_matrix(adj, mode = "directed", weighted = TRUE)
    if (method == "harmonic"){
      temp_d <- 1/igraph::distances(temp_g, mode = mode, algorithm = "dijkstra")
      ## Not consider the distance of a vertex to itself
      diag(temp_d) <- NA
      if (mode == "in"){
        closeness_c_output[,2] <- rowSums(temp_d, na.rm = TRUE)
      }
      if (mode == "out"){
        closeness_c_output[,2] <- rowSums(temp_d, na.rm = TRUE)
      }
    }
    if (method == "standard"){
      temp_d <- igraph::distances(temp_g, mode = mode, algorithm = "dijkstra")
      diag(temp_d) <- NA
      if (mode == "in"){
        closeness_c_output[,2] <- 1/rowSums(temp_d, na.rm = TRUE)
      }
      if (mode == "out"){
        closeness_c_output[,2] <- 1/rowSums(temp_d, na.rm = TRUE)
      }
    }
    return(closeness_c_output)
  }
}

#' Weighted PageRank centrality
#'
#' Compute the weighted PageRank centrality measures of the vertices in a weighted and directed 
#' network represented through its adjacency matrix.
#'
#' @param adj is an adjacency matrix of an weighted and directed network
#' @param gamma is the damping factor; it takes 0.85 (default) if not given.
#' @param theta is a tuning parameter leveraging node degree and strength; theta = 0 does not consider
#' edge weight; theta = 1 (default) fully considers edge weight.
#' @param prior.info vertex-specific prior information for restarting when arriving at a sink. When
#' it is not given (\code{NULL}), a random restart is implemented.
#' 
#' @return a list of vertices with corresponding weighted PageRank scores
#'
#' @references
#' \itemize{
#' \item Zhang, P., Wang, T. and Yan, J. (2021+) PageRank centrality and algorithms for 
#' weighted, directed networks with applications to World Input-Output Tables
#' }
#'
#' @note 
#' Function \code{wpr} is an extension of function \code{page_rank} 
#' in package \code{igraph}.
#'
#' @examples
#' ## Generate a network according to the Erd\"{o}s-Renyi model of order 20
#' ## and parameter p = 0.3
#' edge_ER <- rbinom(400,1,0.3)
#' weight_ER <- sapply(edge_ER, function(x){x*sample(3,1)})
#' adj_ER <- matrix(weight_ER,20,20)
#' mywpr <- wpr(adj_ER, gamma = 0.85, theta = 0.75)
#' system.time(mywpr)
#' 
#' @export

wpr <- function(adj, gamma = 0.85, theta = 1, prior.info){
  
  ## regularity conditions
  if (dim(adj)[1]!=dim(adj)[2]){
    stop("The adjacency matrix is not a square matrix!")
  }
  if ((gamma < 0) | (gamma > 1)){
    stop("The damping factor is not between 0 and 1!")
  }
  if ((theta < 0) | (theta > 1)){
    stop("The tuning parameter is not between 0 and 1!")
  } 
  if (missing(prior.info)){
    prior.info = rep(1/dim(adj)[1],dim(adj)[1])
    warning("No prior information is given; A uniform prior is in use!")
  }
  if (length(prior.info) != dim(adj)[1]){
    stop("The dimension of the prior information is incorrect!")
  }
  if ((sum(prior.info) == 0) | any(prior.info < 0)){
    stop("The prior information is invalid!")
  }
  if (sum(prior.info) != 1){
    prior.info <- prior.info/sum(prior.info)
    warning("The prior information is not normalized!")
  }
  
  ## get the unweighted adjacency matrix
  unweight.adj <- adj
  unweight.adj[unweight.adj > 0] <- 1
  
  ## construct M and M.star matrix
  n <- dim(adj)[1]
  sink.node <- which(rowSums(adj) == 0)
  M <- theta*t(adj/rowSums(adj)) + (1 - theta)*t(unweight.adj/(rowSums(unweight.adj)))
  M[, sink.node] <- prior.info
  B <- matrix(rep(prior.info, n), nrow = n, ncol = n)
  M.star <- gamma*M + (1 - gamma)*B
  
  ## rARPACK cannot solve solve matrices of 2-by-2
  if (dim(adj)[1] == 2){
    eig_sol <- eigen(M.star)
    eigen_v <- eig_sol$vectors[,1]
    eigen_vstd <- abs(eigen_v)/sum(abs(eigen_v))
    name_v <- c(1:n)
    myres <- cbind(name_v, eigen_vstd)
    colnames(myres) <- c("vertex","WPR")
    return(myres)
  }
  
  ## use rARPACK to solve large-scale matrix
  if (dim(adj)[1] > 2){
    eig_sol <-rARPACK::eigs(M.star, k = 1, which = "LM", mattype = "matrix")
    eigen_v <- Re(eig_sol$vectors)
    eigen_vstd <- abs(eigen_v) / sum(abs(eigen_v))
    name_v <- c(1:n)
    myres <- cbind(name_v, eigen_vstd)
    colnames(myres) <- c("vertex","WPR")
    return(myres)
  }
}

