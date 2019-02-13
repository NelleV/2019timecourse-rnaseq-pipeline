library(ClusterR)

#' Performs splines clustering using K-means
#'
#' @param data
#' @param meta
#' @param n_clusters int optional, default: 10
#' @param basis basis, optional, default: NULL 
#'		to use for the splines fitting
#' @param init	["kmeans++", "random", "optimal_init"]
#' @param random_seed int, optional, default: NULL
#' @param degrees_of_freedom int, optional, default: 4
#' @export
splines_kmeans = function(data, meta, n_clusters=10,
			  basis=NULL,
			  init="kmeans++",
			  n_init=10,
			  max_iter=300,
			  random_seed=NULL,
			  fit_splines=TRUE,
			  rescale=TRUE,
			  degrees_of_freedom=4){
    meta = check_meta(meta)
    check_data_meta(data, meta)
    if(is.null(basis)){
        full_model = ~Group:ns(Time, df=df) + Group + 0
	X = model.matrix(full_model, data=meta)
    }else{
	X = basis
    }

    if(fit_splines){
        fitted_data = fit_predict_splines(data, X)
    }else{
	fitted_data = data
    }

    if(rescale){
        fitted_data = rescale_values(fitted_data, meta)
    }

    # Set the random seed if it is null.
    if(is.null(random_seed)){
	set.seed()
	random_seed = .Random.seed[1]
    }
    kmeans_clusters = ClusterR::KMeans_rcpp(
	fitted_data, n_clusters, num_init=n_init, max_iters=max_iter,
	seed=random_seed, initializer=init)
    return(kmeans_clusters)
}
