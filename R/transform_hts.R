#' Summing matrix to Nodes
#'
#' Transform summing matrix to the nodes structure of \code{hts}.
#' @param S summing matrix
#'
#' @return The nodes structure of \code{hts}.
#' @export
#'
#' @examples
#'
#' hierarchy <- rbind(rep(1, 9), c(rep(1, 3), rep(0, 6)), c(rep(0, 3), rep(1, 3), rep(0, 3)), c(rep(0, 6), rep(1, 3)), diag(1, 9))
#' nodes <- Hierarchy2Nodes(hierarchy)
#'
Hierarchy2Nodes <- function(S) {

  size.m <- nrow(S); size.n <- ncol(S)

  indicators.level <- cumsum(rowSums(S)) %% size.n; ends.level <- (1:size.m)[!indicators.level]

  nodes.all <- list(); number.level <- length(ends.level)
  for (i in 2:number.level) {

    if (i == 2)
      S.parent <- t(rbind(S[1:ends.level[i - 1], ]))
    else
      S.parent <- t(S[(ends.level[i - 2] + 1):ends.level[i - 1], ])

    S.child <- t(S[(ends.level[i - 1] + 1):ends.level[i], ])

    nodes.j <- NULL; count.tmp <- 0; S.tmp <- rep(0, size.n)
    for (j in 1:ncol(S.child)) {

      count.tmp <- count.tmp + 1; S.tmp <- S.tmp + S.child[, j]

      test.new <- sum(apply(S.parent == S.tmp, 2, prod))
      if (test.new) {
        nodes.j <- c(nodes.j, count.tmp)
        count.tmp <- 0; S.tmp <- rep(0, size.n)
      }
    }

    nodes.all <- c(nodes.all, list(nodes.j))
  }

  return(nodes.all)
}



#' Nodes to Summing matrix
#'
#' Transform the nodes structure of \code{hts} to summing matrix
#' @param nodes The nodes structure of \code{hts}.
#'
#' @return summing matrix S
#' @export
#'
#' @examples
#' nodes <- list(3,c(3,3,3))
#' S <- Nodes2Hierarchy(nodes)
#'
Nodes2Hierarchy <- function(nodes) {

  number.level <- length(nodes); size.n <- sum(nodes[[number.level]])

  S <- S.child <- diag(1, size.n)
  for (i in number.level:1) {

    S.parent <- NULL
    for (j in nodes[[i]]) {
      S.parent <- rbind(S.parent, colSums(rbind(S.child[1:j, ])))
      S.child <- S.child[-(1:j), ]
    }

    S <- rbind(S.parent, S); S.child <- S.parent
  }

  return(S)
}
