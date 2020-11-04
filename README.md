# A pipeline to analyse transcriptomic time-course data

An automatically compiled version of the manuscript's Rmd with Travis-CI can
be found here: https://nellev.github.io/2019timecourse-rnaseq-pipeline/

1. Install the dependencies and the package.

    - Install all the dependances at once using the following command:
	`Rscript scripts/install.R`
    - Install timecoursedata via the install_github commadd
	install_github("NelleV/timecoursedata", dependencies=FALSE)
    - Install moanin via the install_github command
	install_github("NelleV/moanin", dependencies=FALSE))


2. Compile the manuscript:

    cd scripts
    make

  Or, directly in R:

```r
    library(rmarkdown)
    render("manuscript.Rmd", output_file="reports/manuscript.pdf")
```

3. Other elements: Note that some elements, such as the bootstrapped k-means,
   the consensus plots, and the pathway and GO term enrichment are
   pre-computed. To run everything from scratch:

    - Run the clustering for k between 2 and 20, and bootstrapped random seed
      between 1 and 30:

      * On a cluster, using the script
	`scripts/cluster_scripts/stability_kmeans.sh`
      * On a local machine, with:
        
	```bash
	for SEED in {1..30};
	do
	    for N_CLUSTERS in {2..20};
	    do
		Rscript run_clustering.R ${N_CLUSTERS} ${SEED}
	    done;
	done;
	```

    - Merge all results in a single file using the command:
	```bash
	Rscript concatenate_clustering_stability_results.R
	```
    - Compute the consensus plots:
	```bash
	Rscript run_consensus_plots.R
	```
    - Compute the pathway enrichment analysis for all clusters
	```bash
	Rscript run_pathways.R
	```
    - Compute the GO term enrichment analysis for all clusters
	```bash
	Rscript run_go_terms.R
	```


##  SUPPLEMENTARY: PThe normalization of the micro-array data


You can find an Rmd file containing detailed steps for the normalization in
`scripts/miscellaneous/normalization.Rmd`
