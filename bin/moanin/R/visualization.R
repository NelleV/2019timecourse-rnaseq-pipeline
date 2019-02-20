library(viridis)

#' Plotting centroids
#'
#' @param centroids matrix (k, t) containing the centroids
#' @param meta	data.frame (t, n) containing the metadata.
#' @param colors vector, optional, default NULL
#'		vector of colors
#'
#' @export
plot_centroids = function(centroids, meta, colors=NULL){
    n_centroids = dim(centroids)[1]
    if(n_centroids <= 3){
        par(mfrow=c(n_centroids, 1),
            mar=c(0.1, 0.1, 0.1, 0.1))
    }else if(n_centroids <= 6){
        n_col = ceiling(n_centroids / 2)
        par(mfrow=c(n_col, 2),
            mar=c(0.1, 0.1, 0.1, 0.1))
    }else if(n_centroids <= 12){
        ncol = ceiling(n_centroids / 3)
        par(mfrow=c(ncol, 3),
            mar=c(0.1, 0.1, 0.1, 0.1))
    }else{
        nrow = round(n_centroids ** 0.5)
        ncol = ceiling(n_centroids / nrow)
        par(mfrow=c(ncol, nrow),
            mar=c(0.1, 0.1, 0.1, 0.1))
    }

    for(i in 1:n_centroids){
        plot_centroid_individual(centroids[i, ], meta, colors=colors)
    }
}


plot_centroid_individual = function(centroid, meta, colors=NULL){
    groups = levels(meta$Group)

    xrange = range(meta$Time)
    yrange = range(centroid)
    
    plot(xrange, yrange, type="n")
    if(is.null(colors)){
        colors = viridis(length(groups))
    }
    
    for(i in 1:length(groups)){
        group = groups[i]
        color = colors[i]
        
        mask = meta$Group == group
        time = meta$Time[mask]
        indx = order(time)
        lines(time[indx], centroid[mask][indx], type="b", col=color, pch=16,
              lwd=2)
    }         
} 
