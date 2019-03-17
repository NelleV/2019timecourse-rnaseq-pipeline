# Splines utility

library(splines)
library(MASS)
library(stats)


#' Fit splines
#'
#' @param data the data
#' @param splines_model splines_model 
#' @param weights weigts
#'
#' @return beta coefficients
#'
#' @export
fit_splines = function(data, splines_model, weights=NULL){
    basis = splines_model$basis
    n = ncol(basis)
    nr = nrow(data)
    basis = splines_model$basis 

    if(!is.null(weights)){
	beta = matrix(nrow=nr, ncol=n)
	for(i in 1:nr){
	    beta[i,] = stats::lm.wfit(basis, data[i,], weights[i,])$coefficients
        }
	row.names(beta) = row.names(data)
    }else{
	beta = t(stats::lm.fit(basis, t(data))$coefficients)
    }
    return(beta)
}

#' Fit and predict splines
#'
#' @param data the data
#' @param splines_model splines_model
#' @param weights weigts
#'
#' @return y_fitted the fitted y values
#'
#' @export
fit_predict_splines = function(data, splines_model, weights=NULL){
    basis = splines_model$basis
    y_fitted = t(stats::lm.fit(basis, t(data))$fitted.values)
    return(y_fitted)
}


rescale_values = function(y, meta, group=NULL){
    if(is.null(group)){
	ymin = row_min(y) 
	y = y - ymin
	ymax = row_max(y)
	# We may have a division by 0 here
	y = y / ymax
    }else{
	factors_to_consider = levels(unlist(meta[group]))
	for(factor in factors_to_consider){
	    mask = meta["Group"] == factor
	    ymin = row_min(y[, mask]) 
	    y[, mask] = y[, mask] - ymin
	    ymax = row_max(y[, mask])
	    # We may have a division by 0 here
	    y[, mask] = y[, mask] / ymax
	}
    }
    return(y)
}


# XXX It's wierd that this does not exists in R…
# it probably exists but under another name?
row_max = function(X){
  return(apply(X, 1, max))
}

row_min = function(X){
  return(apply(X, 1, min))
}

row_mean = function(X){
    return(apply(X, 1, mean))
}

row_sum = function(X){
    return(apply(X, 1, sum))
}

row_argmin = function(X){
    return(apply(X, 1, which.min))
}


# Worst name ever
align_data_onto_centroid = function(data, centroid){
    n_samples = dim(data)[2]
    n_genes = dim(data)[1]
    if(n_samples != length(centroid)){
	stop("align_data_onto_centroid: problem in dimensions")
    }
    scaling_factors = (
	rowSums(rep(centroid - mean(centroid), each=n_genes) * data) /
	rowSums((data - rowMeans(data)) * data))
    scaling_factors[scaling_factors < 0] = 0
    shift_factors = rowMeans(
	rep(centroid, each=n_genes) - rep(scaling_factors, times=n_samples) * data)

    data_fitted = rep(scaling_factors, times=n_samples) * data
    data_fitted = (data_fitted + rep(shift_factors, times=n_samples))
    return(data_fitted)
}


score_genes_centroid = function(data, centroid){
    n_genes = dim(data)[1]
    data_fitted = align_data_onto_centroid(data, centroid)

    scores = row_sum((data_fitted - rep(centroid, times=n_genes))**2)

    scores = scores / max(scores)
    return(scores)
}


#' Fisher's method to combine pvalues
#'
#' Combines all p-value per rows.
#' 
#' @param pvalues pvalues
#'
#' @keywords internal
fisher_method = function(pvalues){
    # TODO Add a check that all pvalues are "valid"
    keep = (pvalues >= 0) & (pvalues <= 1)
    pvalues[pvalues == 0] = 1e-285
    
    lnp = log(pvalues)
    chisq = (-2) * row_sum(lnp)
    df = 2 * length(lnp)
    fisher_pval = stats::pchisq(chisq, df, lower.tail=FALSE)
    return(fisher_pval)
}
