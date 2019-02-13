

#' Check that the metadata provided is what we expect
#'
#' This method will raise errors if the metadata provided is not as expected.
#'
#' @param meta metadata
#' @return meta returns the metadata with additional columns if necessary.
check_meta = function(meta){
    metadata_column_names = colnames(meta)
    if(!("Group" %in% metadata_column_names)){
	stop(
	    "Metadata doesn't contain expected information." +
	    " Group column is missing.")
    }

    if(!("Time" %in% metadata_column_names)){
	stop(
	    "Metadata doesn't contain expected information." +
	    " Group column is missing.")
    }

    # Just create this one.
    if(!("WeeklyGroup" %in% metadata_column_names)){
	meta["WeeklyGroup"] = as.factor(make.names(meta$Group:as.factor(meta$Time)))
    }
    return(meta) 
}


#' Check data and meta
check_data_meta = function(data, meta){
    dim_data = dim(data)
    dim_meta = dim(meta)
    data = as.matrix(data)
    if(dim_meta[1] != dim_data[2]){
	stop(
	    "Data and metadata are inconsistent. Data is of shape (Xx"+
	    "Metadata is of shape XX")
    }

    if(!is.numeric(data)){
	stop("Data should be of type numeric")
    }
    meta = check_meta(meta)
}
