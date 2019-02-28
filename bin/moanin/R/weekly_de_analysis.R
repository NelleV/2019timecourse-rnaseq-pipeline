library("limma")
library("stats")
library("edgeR")


#' Fit weekly differential expression analysis
#'
#' @param counts Gene expression data
#' @param meta Metadata data frame
#' @param contrasts Contrast to use.
#' @param use_voom_weights boolean: whether to use voom weights or not
#'
#' @export
weekly_differential_expression_analysis = function(counts, meta,
						   contrasts,
						   use_voom_weights=TRUE){
    meta = check_meta(meta)

    design = stats::model.matrix(~WeeklyGroup + 0, data=meta)

    cleaned_colnames = gsub("WeeklyGroup", "", colnames(design))
    colnames(design) = cleaned_colnames

    allcontrasts = limma::makeContrasts(
        contrasts=contrasts,
        levels=design)

    if(use_voom_weights){
        y = edgeR::DGEList(counts=counts)
	y = edgeR::calcNormFactors(y, method="upperquartile")
        v = limma::voom(y, design, plot=FALSE)
	v = limma::lmFit(v)
    }else{
	v = limma::lmFit(counts, design)	
    }

    fit = limma::contrasts.fit(v, allcontrasts)
    fit = limma::eBayes(fit)
    contrast_names = colnames(fit$p.value)
    fit$adj.p.value = stats::p.adjust(fit$p.value, method="BH")
    dim(fit$adj.p.value) = dim(fit$p.value)
    colnames(fit$adj.p.value) = contrast_names

    combine_results = function(ii, fit2){
	contrast_formula = contrasts[ii]
	de_analysis = data.frame(row.names=row.names(counts))

        base_colname = gsub(" ", "", contrast_formula, fixed=TRUE)
	colname_pval = paste(base_colname, "-pval", sep="")
	colname_qval = paste(base_colname, "-qval", sep="")
	colname_lfc = paste(base_colname, "-lfc", sep="")

        tt = limma::topTable(
            fit2, coef=ii, number=length(rownames(fit2$coef)),
            p.value=1, adjust.method="none",
            genelist=rownames(fit2$coef))
	de_analysis[colname_pval] = fit2$p.value[, contrast_formula]
	de_analysis[colname_qval] = fit2$adj.p.value[,  contrast_formula]
	de_analysis[colname_lfc] = tt$logFC
	return(de_analysis)
    }

    all_results = do.call("cbind",
			  lapply(1:length(contrast_names),
				 combine_results, fit2=fit))
    return(all_results)
}
