
# This should override the Biobase function
rowMax = function(data){
    return(as.vector(lapply(data, function(x) max(x))))
}

rowMin = function(data){
    return(as.vector(lapply(data, function(x) min(x))))
}

