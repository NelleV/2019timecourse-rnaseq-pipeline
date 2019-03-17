library(ClusterR)
library(stats)
library(splines)

#' Performs splines clustering using K-means
#'
#' @param data matrix containg the data. Data needs to be in log scale.
#' @param splines_model splines_model 
#' @param n_clusters int optional, default: 10
#' @param init	["kmeans++", "random", "optimal_init"]
#' @param n_init int, optional, default: 10
#'	    Number of initialization to perform.
#' @param max_iter  int, optional, default: 300
#'	Maximum number of iteration to perform	
#' @param random_seed int, optional, default: NULL
#' @param fit_splines	boolean, optional, default: TRUE
#'	Whether to fit splines or not.
#' @param rescale   boolean, optional, default: TRUE
#'	Whether to rescale the data or not.
#' @export
splines_kmeans = function(data, splines_model, n_clusters=10,
			  init="kmeans++",
			  n_init=10,
			  max_iter=300,
			  random_seed=NULL,
			  fit_splines=TRUE,
			  rescale=TRUE){
    meta = splines_model$meta
    basis = splines_model$basis
    check_data_meta(data, meta)

    if(fit_splines){
        fitted_data = fit_predict_splines(data, splines_model)
    }else{
	fitted_data = data
    }

    if(rescale){
        fitted_data = rescale_values(fitted_data, meta)
    }

    # Set the random seed if it is null.
    if(is.null(random_seed)){
	set.seed(NULL)
	random_seed = .Random.seed[1]
    }
    kmeans_clusters = ClusterR::KMeans_rcpp(
	fitted_data, n_clusters, num_init=n_init, max_iters=max_iter,
	seed=random_seed, initializer=init)

    return(kmeans_clusters)
}
