library(topGO)


# ontology: BP, CC, NF
find_enriched_go_terms = function(scores, geneID2GO,
				  ontology="BP", 
				  weighted=FALSE,
				  node_size=10){
    if(!(ontology %in% c("BP", "CC", "NF"))){
	error_message = paste(
	    "moanin::find_enriched_go_terms: Ontology should be 'BP', 'CC',",
	    "or 'NF'. Ontology provided is",
	    ontology)
	stop(error_message)
    }


    getTopDiffGenes = function(data, cutOff=NULL){
      return(data < 0.5)
    }

    GOdata = new("topGOdata", ontology=ontology, allGenes=scores, nodeSize=node_size,
		 geneSel=getTopDiffGenes,
		 annot=annFUN.gene2GO, gene2GO=geneID2GO)

    if(weighted){
        resultFisher = runTest(GOdata, algorithm="classic", statistic="fisher")
    }else{
        resultFisher = runTest(GOdata, algorithm="weight", statistic="fisher")
    }
    n_nodes = length(resultFisher@score)

    allRes = GenTable(
	GOdata, 
	resultFisher=resultFisher,
	orderBy="resultFisher", ranksOf="resultFisher",
	topNodes=n_nodes)
        
    # P-value correct
    allRes[, "resultFisher_padj"] = stats::p.adjust(allRes$resultFisher, method="BH")
    wh = which(allRes[, "resultFisher_padj"] <= 0.05)

    allRes = allRes[wh,]
    wh = which(apply(allRes[, c("Significant", "Expected")],
    	              1, function(x){x["Significant"] > x["Expected"]}))
    allRes = allRes[wh,]

    return(allRes)
}
