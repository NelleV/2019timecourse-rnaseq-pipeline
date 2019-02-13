library("moanin")

context("moanin::splines.R")

test_that("splines:align_data_onto_centroid", {
    set.seed(42)
    n_samples = 20
    n_genes = 5
    centroid = runif(n_samples)

    shift = runif(n_genes)
    scale = 1 + runif(n_genes)

    data = rep(centroid, each=n_genes)
    dim(data) = c(n_genes, n_samples)
    expect_equal(data, align_data_onto_centroid(data, centroid))

    # Ok, now let's make this a bit more complicated.
    scale = rep(scale, times=n_samples)
    dim(scale) = dim(data)
    shift = rep(shift, times=n_samples)
    dim(shift) = dim(data)

    shifted_scaled_data = scale * data + shift
    shifted_scaled_centroid = scale * centroid + shift
    expect_equal(data,
		 align_data_onto_centroid(shifted_scaled_data, centroid))
    expect_error(align_data_onto_centroid(data[, 1:10], centroid))
})

