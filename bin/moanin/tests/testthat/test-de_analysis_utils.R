library("moanin")

context("moanin::de_analysis_utils.R")


test_that("Estimating log fold change smoke tests", {
    # This is just a smoke test.
    data(shoemaker2015)
    data = shoemaker2015$data
    meta = shoemaker2015$meta

    # Reduce the data set
    data = data[1:10, ]
    contrasts = c("C-K")
    methods = moanin:::ALL_LFC_METHODS
    for(method in methods){
        expect_silent(
	    estimate_log_fold_change(
		data, meta, contrasts, method=method))
    }
})

test_that("Estimating log fold change with unknown error", {
    data(shoemaker2015)
    data = shoemaker2015$data
    meta = shoemaker2015$meta

    # Reduce the data set
    data = data[1:10, ]
    contrasts = c("C-K")
    expect_error(estimate_log_fold_change(data, meta, contrast, method="hahaha"))

})
