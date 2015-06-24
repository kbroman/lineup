context("corbetw2mat")

test_that("corbetw2mat works with paired=TRUE", {

    n_ind <- 20
    n_col <- 3

    x <- matrix(rnorm(n_ind*n_col), ncol=n_col)
    y <- x + matrix(rnorm(n_ind*n_col, 0, 0.5), ncol=n_col)
    colnames(x) <- colnames(y) <- 1:n_col
    rownames(x) <- rownames(y) <- 1:n_ind

    result <- corbetw2mat(x, y, what="paired")

    expected <- NULL
    for(i in 1:ncol(x))
        expected[i] <- cor(x[,i], y[,i])
    names(expected) <- 1:n_col

    expect_equal(result, expected)

})
