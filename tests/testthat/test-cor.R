context("Correlation")

test_that("sort_vector_abs", {
  #x <- c(-2, -10, 3, 11)
  #sort_vector_abs(x)
  
  set.seed(1)
  for (i in 1:10) {
    x <- runif(100, min = -10, max = 10)
    x_ord <- order(abs(x), decreasing = TRUE)
    expect_equal(x[x_ord], sort_vector_abs(x))
  }
  
  if (FALSE) {
    set.seed(1)
    x <- runif(10000, min = -10, max = 10)
    microbenchmark::microbenchmark(
      r = {
        x_ord <- order(abs(x), decreasing = TRUE)
        y <- x[x_ord]
        force(y)
      },
      rcpp = sort_vector_abs(x),
      times = 100
    )
  }
})

data(danedb)

dmat <- apply(as.matrix(danedb), 2, as.numeric)
res <- pearson_correlation_absolute_sparse_all(dmat, progress = FALSE)
cormat <- cor(dmat)

test_that("pearson_correlation_absolute_sparse_all", {
  expect_true(res$error == FALSE)
  
  # Correlation calculations
  for (i in seq_len(nrow(res$idx))) {
    c1 <- cor(dmat[, res$idx[i, 1]], dmat[, res$idx[i, 2]])
    c2 <- res$corr[i]
    expect_equal(c1, c2)
  }
  
  # Correlation ordering
  cors <- cormat[upper.tri(cormat)]
  cors_ord <- order(abs(cors), decreasing = TRUE)
  expect_equal(cors[cors_ord], res$corr)
})
