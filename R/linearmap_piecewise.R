#  nucdataBaynet - Nuclear Data Evaluation Using a Bayesian Network 
#  Copyright (C) 2019  Georg Schnabel
#  
#  nucdataBaynet is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  nucdataBaynet is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>


#' Create Piecewise Linear Map
#'
#' @param expEn TODO
#' @param opts TODO
#' @param ret.mat TODO
#'
#' @return TODO
#' @import data.table Matrix methods
#' @export
#'
pwMap <- function(expEn, opts = NULL, ret.mat = TRUE) {

  stopifnot(is.list(opts))
  order <- if (is.null(opts$order)) 1 else opts$order
  inv <- if (is.null(opts$inv)) FALSE else opts$inv
  outsideZero <- if (is.null(opts$outsideZero)) FALSE else opts$outsideZero
  srcEn <- opts$ens

  stopifnot(order %in% c(1,2,3))
  stopifnot(is.numeric(srcEn))
  stopifnot(!is.unsorted(srcEn))

  rowIdx <- seq_along(expEn)
  colIdx1 <- findInterval(expEn, srcEn, rightmost.closed = TRUE)
  colIdx2 <- colIdx1 + 1
  if (!isTRUE(outsideZero)) {
    if (!(all(colIdx1 > 0 & colIdx1 < length(srcEn))))
      stop("Some energies in expEn are outside the grid given in opts$ens")
  }

  sel <- colIdx1 > 0 & colIdx2 <= length(srcEn)
  rowIdx <- rowIdx[sel]
  colIdx1 <- colIdx1[sel]
  colIdx2 <- colIdx2[sel]

  len <- srcEn[colIdx2] - srcEn[colIdx1]
  x1 <- (srcEn[colIdx2] - expEn[sel]) / len
  x2 <- (expEn[sel] - srcEn[colIdx1]) / len

  S1 <- sparseMatrix(i = rep(rowIdx, 2), j = c(colIdx1, colIdx2),
                     x = c(x1, x2), dims = c(length(expEn), length(srcEn)))
  S2 <- as(tril(matrix(1, length(srcEn), length(srcEn))), "sparseMatrix")

  tmpS3 <- tril(matrix(1, length(srcEn)-1, length(srcEn)-1))
  S3 <- Diagonal(x = rep(1, length(srcEn)))
  S3[-1,-1] <- tmpS3

  if (order == 1) S <- drop0(S1)
  else if (order == 2) S <- drop0(S1 %*% S2)
  else if (order == 3) S <- drop0(S1 %*% S2 %*% S3)
  else stop("invalid order specification")
  if (!inv) S else solve(S)
  S <- as(S, "sparseMatrix")
  if (ret.mat) S else {
    tmp <- summary(S)
    data.table(IDX1 = tmp$i, IDX2 = tmp$j, X = tmp$x)
  }
}
