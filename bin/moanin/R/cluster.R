library(ClusterR)
library(stats)
library(splines)

#' Performs splines clustering using K-means
#'
#' @param data matrix containg the data. Data needs to be in log scale.
#' @param meta data.frame containing the metadata. Metadata needs to contain
#'	       column 'Group' and 'Time'
#' @param n_clusters int optional, default: 10
#' @param basis basis, optional, default: NULL 
#'		to use for the splines fitting
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
#' @param degrees_of_freedom int, optional, default: 4
#'	Number of degrees of freedom used in the clustering. Only used if
#'	fit_splines is TRUE
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
        full_model = ~Group:splines::ns(Time, df=degrees_of_freedom) + Group + 0
	X = stats::model.matrix(full_model, data=meta)
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
	set.seed(NULL)
	random_seed = .Random.seed[1]
    }
    kmeans_clusters = ClusterR::KMeans_rcpp(
	fitted_data, n_clusters, num_init=n_init, max_iters=max_iter,
	seed=random_seed, initializer=init)
    return(kmeans_clusters)
}
