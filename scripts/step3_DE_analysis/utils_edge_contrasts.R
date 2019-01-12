# This file contains EDGE code adapted to work with limma contrasts

library(edge)
library(splines)
library(MASS)
source("utils_splines.R")

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
  stat = 0
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
  return(stat)
}

compute_pvalue = function(X, y, beta, beta_null, ng_labels, statistics="lrt",
			  df2=NULL, weights=NULL, mask=NULL, developmental=FALSE){

    fitFull = beta %*% t(X)
    if(developmental){
      fitNull = row_mean(y[, mask])
    }else{
      fitNull = beta_null %*% t(X)
    }
  
    if(!is.null(mask)){
      fitFull = fitFull[, mask]
      if(!developmental){
        fitNull = fitNull[, mask]
      }
      y = y[, mask]
      ng_labels = ng_labels[mask]
    }

    if(!is.null(weights)){
      resNull = weights^(1/2) * (y - fitNull)
      resFull = weights^(1/2) * (y - fitFull)
    }else{
      resNull = y - fitNull
      resFull = y - fitFull
    }

    ng = nlevels(ng_labels)
    df = ncol(X) / ng

    if(statistics == "ftest"){
      stat = edge:::lrtStat(resNull, resFull)
      if(is.null(df2)){
        df2 = nrow(X) - ncol(X)
      }
      df1 = df
      pval = 1 - pf(stat * df2 / df1, df1=df1, df2=df2)

    }else{
      lstat = lrtStat(resNull, resFull, ng_labels=ng_labels)
      pval = 1 - pchisq(lstat, df=df)
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
#' @param contrast Contrast using \code{makeContrasts} from \code{limma}.
#' @param center: boolean
#' @param weights Weights matrix for the fitting.
#' @param df Degrees of freedom to for the splines.
#' @param basis the basis matrix. If provided, ignore some of the parameters
#' @param mask boolean mask of which points to use for the tests. This is
#'	  useful when you want to use all data points for fitting the splines,
#'	  but not all the data points for the tests.
edgeWithContrasts = function(data, meta, contrasts=NULL, center=FALSE, weights=NULL, df=4,
			     basis=NULL, mask=NULL, developmental=FALSE){
  ng = nlevels(meta$Group)
  ng_labels = meta$Group

  if(is.null(contrasts) & !developmental){
    # FIXME need better error message
    stop("Needs either contrasts or developmental")
  }

  contrasts_coef = c(contrasts)
  meta = droplevels(meta)

  if(!is.null(mask)){
    if(length(mask) != dim(data)[2]){
      stop("The mask provided doesn't have the correct length.")
    }
  }

  if(length(contrasts_coef) != ng){
    stop("The contrast coef vector should be of the same size as the number of groups")
  }
  

  if(is.null(basis)){
    full_model = ~Group:ns(Time, df=df) + Group + 0
    X = model.matrix(full_model, data=meta)
  }else{
    X = basis
  }
  y = data

  if(center){
    y = center_data(y)
    X = t(center_data(t(X)))
  }


  beta = fit_splines(y, X, weights=weights)
  if(developmental){
    # For developmental, only consider the group from the contrast
    beta_full = compute_beta_null(X, beta, contrasts_coef)
    beta_null = NULL
  }else{
    beta_null = compute_beta_null(X, beta, contrasts_coef)
  }

  pval = compute_pvalue(X, y, beta, beta_null, ng_labels, weights=weights,
		        mask=mask,
		        developmental=developmental)
  pval_BH = p.adjust(pval, method="BH")
  pval_ftest = compute_pvalue(
      X, y, beta, beta_null, ng_labels,
      statistics="ftest",
      mask=mask,
      developmental=developmental)
  pval_ftest_BH = p.adjust(pval_ftest, method="BH")

  fit = NULL
  fit$pval_lrt = pval
  fit$pval = pval_ftest
  fit$qval = pval_ftest_BH
  fit$qval_lrt = pval_BH
  fit$beta = beta
  fit$beta_null = beta_null
  if(is.null(basis)){
    fit$full_model = full_model
  }else{
    fit$full_model = NULL
  }
  fit$basis = X
  fit$weights = weights
  fit$center = center # Keep track of that for when predicting
  fit$contrasts = contrasts
  return(fit)
}


predict_ = function(fit, meta, time, null_model=FALSE){
  meta_predict = expand.grid(
      Group=levels(meta$Group),
      Time=time)

  if(!is.null(fit$full_model)){
    X_pred = model.matrix(fit$full_model, data=meta_predict)
  }else{
    X_pred = fit$basis
    meta_predict = meta  
  }
  
  if(fit$center){
    X_pred = t(center_data(t(X_pred)))
  }
  
  if(null_model){
    y_fitted = fit$beta_null %*% t(X_pred)
  }else{
    y_fitted = fit$beta %*% t(X_pred)
  }

  colnames(y_fitted) = make.names(meta_predict$Group:as.factor(meta_predict$Time))
  return(y_fitted)
}


# FIXME this should not be done this way.
compute_lfc = function(fit, meta, method="mean"){
  
  if(!(method %in% c("mean", "max"))){
    stop("Unknown method to compute log fold change")
  }

  meta = droplevels(meta)
  time_min = min(as.numeric(meta$Time))
  time_max = max(as.numeric(meta$Time))
  nobs = 100
  time = seq(from=time_min, to=time_max, length.out=nobs)

  if(is.null(fit$full_model)){
    meta_predict = meta
    time = meta$Time
    nobs = length(unique(meta$Time))
  }else{
    meta_predict = expand.grid(
        Group=levels(meta$Group),
	Time=time)
  }
  y_fitted = predict_(fit, meta, time)
  y_fitted = t(t(y_fitted) * fit$contrasts[meta_predict$Group])
  y_fitted = summarise(y_fitted, meta_predict$Group:as.factor(meta_predict$Time))
  
  for(i in seq(1:(dim(y_fitted)[2]/nobs))){
    if(i == 1){
	foldchange = y_fitted[, ((i-1)*nobs+1):(i*nobs)]
    }else{
	foldchange = foldchange + y_fitted[, ((i-1)*nobs+1):(i*nobs)]
    } 
  }

  if(any(is.na(foldchange))){
    warning("NAs in the foldchange")
  }

  foldchange[is.na(foldchange)] = 0
  foldchange_ = rowMeans(abs(foldchange))
  fit$lfc_mean = log2(foldchange_+1)
  foldchange = row_max(abs(foldchange))
  fit$lfc_max = log2(foldchange+1)

  if(method == "mean"){
    fit$lfc = fit$lfc_mean
  }else{
    fit$lfc = fit$lfc_max
  }
  return(fit)
}

