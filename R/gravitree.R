#' @export
gravitree <- function(x, m, k = NULL, sample = 1, random = TRUE, na_rm = FALSE) {
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }

  nr <- nrow(x)

  if (nr <= 3) {
    stop("'x' must contain at least 4 observations.")
  }

  if (!is.matrix(m)) {
    m <- as.matrix(m)
  }

  if (nr != nrow(m)) {
    stop("The length of 'm' does not equal the number of rows of 'x'.")
  }

  if (na_rm) {
    na_ix <- apply(x, 1, function(r) any(is.na(r))) | is.na(m)
    x <- x[!na_ix, , drop = FALSE]
    m <- m[!na_ix, , drop = FALSE]
  }

  if (sample < 1) {
    if (random) {
      ix <- sample(1:nr, nr * sample, prob = m)
    } else {
      ix <- order(m, decreasing = TRUE)[1:(nr * sample)]
    }

    notix <- which(!(1:nr %in% ix))
    xx <- x[ix, , drop = FALSE]
    mm <- m[ix, , drop = FALSE]
    g <- Rfast::Tcrossprod(mm, mm) / (Rfast::Dist(xx)^ncol(xx))
    d <- .as_dist(log(1 / g))
  } else {
    g <- Rfast::Tcrossprod(m, m) / (Rfast::Dist(x)^ncol(x))
    d <- .as_dist(log(1 / g))
  }

  tree <- fastcluster::hclust(d, method = "ward.D2")

  if (is.null(k)) {
    k <- length(tree$height) - which.max(diff(tree$height)) + 1
  }

  cluster <- stats::cutree(tree, k)

  if (sample < 1) {
    centers <- lapply(1:k, function(i) {
      ix <- cluster == i
      .wcov(xx[ix, , drop = FALSE], mm[ix])
    })
  } else {
    centers <- lapply(1:k, function(i) {
      ix <- cluster == i
      .wcov(x[ix, , drop = FALSE], m[ix])
    })
  }

  cluster <- Rfast::rowMins(sapply(1:k, function(i) {
    .Mahalanobis(x, centers[[i]]$center, centers[[i]]$cov)
  }))

  list(
    tree = tree,
    cluster = cluster,
    centers = centers
  )
}