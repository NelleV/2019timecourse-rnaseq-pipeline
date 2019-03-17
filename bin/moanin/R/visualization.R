library(viridis)
library(graphics)

#' Plotting centroids
#'
#' @param centroids matrix (k, t) containing the centroids
#' @param splines_model splines_model
#' @param meta	data.frame (t, n) containing the metadata.
#' @param colors vector, optional, default NULL
#'		vector of colors
#'
#' @export
plot_centroids = function(centroids, splines_model, colors=NULL){
    n_centroids = dim(centroids)[1]
    if(n_centroids <= 3){
        graphics::par(
	    mfrow=c(n_centroids, 1),
            mar=c(0.1, 0.1, 0.1, 0.1))
    }else if(n_centroids <= 6){
        n_col = ceiling(n_centroids / 2)
        graphics::par(
	    mfrow=c(n_col, 2),
            mar=c(0.1, 0.1, 0.1, 0.1))
    }else if(n_centroids <= 12){
        ncol = ceiling(n_centroids / 3)
        graphics::par(
	    mfrow=c(ncol, 3),
            mar=c(0.1, 0.1, 0.1, 0.1))
    }else{
        nrow = round(n_centroids ** 0.5)
        ncol = ceiling(n_centroids / nrow)
        graphics::par(
	    mfrow=c(ncol, nrow),
            mar=c(0.1, 0.1, 0.1, 0.1))
    }

    for(i in 1:n_centroids){
        plot_centroid_individual(centroids[i, ], splines_model, colors=colors)
    }
}


plot_centroid_individual = function(centroid, splines_model, colors=NULL){
    meta = splines_model$meta
    groups = levels(meta$Group)

    xrange = range(meta$Timepoint)
    yrange = range(centroid)
    
    graphics::plot(xrange, yrange, type="n")
    if(is.null(colors)){
        colors = viridis::viridis(length(groups))
    }
    
    # scatter points for values
    for(i in 1:length(groups)){
        group = groups[i]
        color = colors[i]
        
        mask = meta$Group == group
        time = meta$Timepoint[mask]
        indx = order(time)
        graphics::lines(time[indx], centroid[mask][indx], type="b",
		        col=color, pch=16,
			lwd=0)
    }
} 


plot_gene_splines = function(data, meta, gene_name, colors=NULL){
    # First, select the gene:
    if(!(gene_name %in% row.names(data))){
	msg = paste("moanin::plot_gene_splines: The gene_name provided '",
		    gene_name, "' is not in the data", sep="")
    }
    gene_data = data[gene_name, ]

    # All of this should be extracted from the amazing model object that we
    # don't have implemented yet.
    degrees_of_freedom = 6
    model = ~Group:splines::ns(Timepoint, degrees_of_freedom) + Group + 0
    
}
