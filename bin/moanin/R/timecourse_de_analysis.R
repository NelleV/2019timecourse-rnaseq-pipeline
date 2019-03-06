# This file contains EDGE code adapted to work with limma contrasts

center_data = function(y, ng_labels){
  for(g in levels(ng_labels)){
    whKeep = which(ng_labels == g)
    sub_mean = rowMeans(y[, whKeep])
    y[, whKeep] = y[, whKeep] - sub_mean
  }
  return(y)
}

compute_beta_null = function(X, beta, contrasts_coef){
  ng = length(contrasts_coef)
  df = ncol(X) / ng
  contrasts_coef_ = rep(contrasts_coef, times=df)

  # Reshape b so that each row corresponds to a group and drop the intercept
  b_ = array(beta, dim=c(dim(beta)[1], ng, df))
  # First start by constructing the matrix T
  t_ = apply(b_, 1, function(x){colSums(x*contrasts_coef)})

  # The observations can not be assumed to be balanced...
  # We need to get rid of the intercept for this part
  # FIXME don't invert this matrix…
  part_K = MASS::ginv(t(X) %*% X)
  K = part_K * contrasts_coef_**2

  # We now need to sum all elements associated to the same pairs of splines.
  # Which is, in a particular case every four elements in both directions.

  K = sapply(1:df, function(jg) rowSums(K[, (jg-1)*ng + 1:ng]))
  K = sapply(1:df, function(jg) colSums(K[(jg-1)*ng + 1:ng,]))

  T_ = MASS::ginv(K) %*% t_

  # We got T. Now, let's move on to the rest
  tmp = as.array(rep(as.vector(T_), each=ng), dim=c(1, 1, 1))
  dim(tmp) = c(ng * df, dim(beta)[1])
  C_ = contrasts_coef * part_K %*% tmp

  C_ = t(C_)
  dim(C_) = c(dim(beta)[1], ng, df)
  beta_null = b_ - C_

  # Last step: reshape beta null so that it is of the same shape as beta
  dim(beta_null) = dim(beta)
  return(beta_null)
}


lrtStat = function(resNull, resFull, ng_labels=NULL) {
  # FIXME I'm pretty sure that in the case of contrasts, the degrees of
  # freedom computed here are wrong as they include part of the data that is
  # not used for the test. This needs to be fixed
  stat = 0
  if(is.null(ng_labels)){
    ss0 = rowSums(resNull^2)
    ss1 = rowSums(resFull^2)
    n = ncol(resNull)
    stat = stat + n * (ss0 - ss1)/(ss1)

  }else{
    for(g in levels(ng_labels)){
	whKeep = which(ng_labels == g)
	sub_resNull = resNull[,whKeep]
	sub_resFull = resFull[,whKeep]

	# Somehow the two lines above don't return the same object depending on
	# the dimension of resNull and resFull, so need to distinguish the case
	# where there is only one observation in data.
	if(is.null(dim(sub_resNull))){
	    ss0 = sum(sub_resNull^2)
	    ss1 = sum(sub_resFull^2)
	    n = length(sub_resNull)
	}else{
	    ss0 = rowSums(sub_resNull^2)
	    ss1 = rowSums(sub_resFull^2)
	    n = ncol(sub_resNull)
	}
	stat = stat + n * (ss0 - ss1)/(ss1)
    }
  }

  return(stat)
}

compute_pvalue = function(X, y, beta, beta_null, ng_labels,
			  n_groups=NULL,
			  n_samples=NULL,
			  degrees_of_freedom=NULL,
			  statistics="lrt",
			  df2=NULL, weights=NULL){

    fitFull = beta %*% t(X)

    fitNull = beta_null %*% t(X)

    if(!is.null(weights)){
      resNull = weights^(1/2) * (y - fitNull)
      resFull = weights^(1/2) * (y - fitFull)
    }else{
      resNull = y - fitNull
      resFull = y - fitFull
    }

    # estimate degrees of freedom.
    if(is.null(n_groups)){
        n_groups = nlevels(ng_labels)
	# FIXME Raise warning
    }

    if(is.null(n_samples)){
	n_samples = ncol(X)
	# FIXME raise warning
    }
    
    df = degrees_of_freedom

    if(statistics == "ftest"){
	stat = lrtStat(resNull, resFull)
	if(is.null(df2)){
	    # FIXME Check this.
	    df2 = n_samples - degrees_of_freedom * n_groups
	}
	df1 = df
	pval = stats::pf(stat * df2 / df1, df1=df1, df2=df2, lower.tail=FALSE)
    }else{
	lstat = lrtStat(resNull, resFull, ng_labels=ng_labels)
	pval = stats::pchisq(lstat, df=degrees_of_freedom, lower.tail=FALSE)
    }
    return(pval)
}

summarise = function(X, ng_levels) {
  X_mean = matrix(,nrow=nrow(X), ncol=nlevels(ng_levels))
  colnames(X_mean) = levels(ng_levels)
  rownames(X_mean) = rownames(X)

  for(g in levels(ng_levels)){
    whKeep = which(ng_levels == g)
    if(length(whKeep) > 1){
      X_mean[, g] = rowMeans(X[, whKeep])
    }else if(length(whKeep) != 0){
      X_mean[, g] = X[, whKeep]
    }
  }
  return(X_mean)
}

#' Run edge with contrasts.
#' 
#' @param data The data matrix.
#' @param meta Meta data which should be consistent with the column names \code{data}.
#' @param contrasts Contrast using \code{makeContrasts} from \code{limma}.
#' @param center boolean, whether to center the data matrix
#' @param weights Weights matrix for the fitting.
#' @param df Degrees of freedom to for the splines.
#' @param basis the basis matrix. If provided, ignore some of the parameters
#' @export
timecourse_differential_expression_analysis = function(data,
						       meta,
						       contrasts,
						       center=FALSE,
						       weights=NULL,
						       df=4,
						       basis=NULL){
    ng = nlevels(meta$Group)
    ng_labels = meta$Group

    meta = check_meta(meta)
    meta = droplevels(meta)
    check_data_meta(data, meta)


    if(is.null(basis)){
	full_model = ~Group:splines::ns(Time, df=df) + Group + 0
	X = stats::model.matrix(full_model, data=meta)
    }else{
	X = basis
    }
    y = data

    if(center){
	y = center_data(y)
	X = t(center_data(t(X)))
    }

    beta = fit_splines(y, X, weights=weights)

    if(dim(contrasts)[1] != ng){
	stop("The contrast coef vector should be of the same size" +
	     " as the number of groups")
    }

    pvalues = data.frame(row.names=row.names(data))
    for(col in 1:ncol(contrasts)){
	contrast = contrasts[, col]

	# Create the name of the column
	contrast_name = colnames(contrasts)[col]
	contrast_name = gsub(" ", "", contrast_name, fixed=TRUE)

	# Get the number of samples used for this particular contrast:
	groups_of_interest = names(contrast)[contrast != 0]
	n_samples_fit = sum(with(meta, Group %in% groups_of_interest))
	n_groups = length(groups_of_interest)
	degrees_of_freedom = dim(X)[2] / ng

	beta_null = compute_beta_null(X, beta, contrast)

	pval = compute_pvalue(X, y, beta, beta_null, ng_labels, weights=weights,
			      n_samples=n_samples_fit,
			      n_groups=n_groups,
			      degrees_of_freedom=degrees_of_freedom)
	pvalues[contrast_name] = pval
    }
    return(pvalues)
}


