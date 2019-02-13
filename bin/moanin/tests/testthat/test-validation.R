library(moanin)

context("Testing validation functions")

test_that("validation:check_meta", {
    data(shoemaker2015)
    data = shoemaker2015$data
    meta = shoemaker2015$meta

    data = data[1:10, ]

    # Running check_data_meta should work fine
    check_data_meta(data, meta)

    expect_error(
	check_data_meta(data[, 2:10], meta))
})
