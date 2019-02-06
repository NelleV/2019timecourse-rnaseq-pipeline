# Splines utility

library(splines)
library(MASS)


#' Fit splines
#'
#' @param y the data
#' @param X the basis
#' @param weights weigts
#'
#' @return beta coefficients
#'
#' @export
fit_splines = function(y, X, weights=NULL){
    n = ncol(X)
    nr = nrow(y)

    
    if(!is.null(weights)){
	beta = matrix(nrow=nr, ncol=n)
	for(i in 1:nr){
	    beta[i,] = lm.wfit(X, y[i,], weights[i,])$coefficients
        }
	row.names(beta) = row.names(data)
    }else{
	# Don't inverse directly the matrix
	beta = t(lm.fit(X, t(y))$coefficients)
    }
    return(beta)
}

#' Fit and predict splines
#'
#' @param y the data
#' @param X the basis
#' @param weights weigts
#'
#' @return y_fitted the fitted y values
#'
#' @export
fit_predict_splines = function(y, X, weights=NULL){
    y_fitted = t(lm.fit(X, t(y))$fitted.values)
    return(y_fitted)
}


rescale_values = function(y, meta, group="Group"){
    factors_to_consider = levels(unlist(meta[group]))
    for(factor in factors_to_consider){
	mask = meta["Group"] == factor
	ymin = row_min(y[, mask]) 
	y[, mask] = y[, mask] - ymin
	ymax = row_max(y[, mask])
	# We may have a division by 0 here
	y[, mask] = y[, mask] / ymax
    }

    return(y)
}


# XXX It's wierd that this does not exists in R…
# it probably exists but under another name?
row_max = function(X){
  return(apply(X, 1, max))
}


# XXX It's wierd that this does not exists in R…
# it probably exists but under another name?
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

align_data_onto_centroid = function(y, centroid){
    scaling_factors = (
	rowSums(centroid - mean(centroid) * y) /
	rowSums(y - row_mean(y) * y))
    scaling_factors[scaling_factors < 0] = 0
    shift_factors = row_mean(centroid - scaling_factors * y)
    return(scaling_factors * y + shift_factors)
}


score_genes_centroid = function(y, centroid){
    y_fitted = align_data_onto_centroid(y, centroid)
    scores = row_sum((y_fitted - centroid)**2)
    scores = scores / max(scores)
    return(scores)
}


#' Fisher's method to combine pvalues
#'
#' Combines all p-value per rows.
#' 
#' @param pvalues pvalues
fisher_method = function(pvalues){
    # TODO Add a check that all pvalues are "valid"
    keep = (pvalues >= 0) & (pvalues <= 1)
    pvalues[pvalues == 0] = 1e-285
    
    lnp = log(pvalues)
    chisq = (-2) * row_sum(lnp)
    df = 2 * length(lnp)
    fisher_pval = pchisq(chisq, df, lower.tail=FALSE)
    return(fisher_pval)
}
