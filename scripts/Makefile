RMD_FILES = $(wildcard *.Rmd)
HTML_FILES = $(patsubst %.Rmd, reports/%.html, $(RMD_FILES))
R_FILES = $(patsubst %.Rmd, tmp/%.R, $(RMD_FILES))

PIPELINE_STEPS = step0_preprocessing \
		 step2_normalization \
		 step3_DE_analysis \
		 step4_clustering

all: html

html: $(HTML_FILES)

complete: html
	make results/all_go_terms.tsv
	rm -f reports/04_downstream_analysis.html
	make html

pipeline:
	for i in $(PIPELINE_STEPS); do $(MAKE) -C $$i all; done


reports/02_differential_expression_analysis.html: 02_differential_expression_analysis.Rmd reports/01_quality_control_normalization.html
	mkdir -p reports
	Rscript .scripts/compile_rmarkdown.R $< $@

reports/03_clustering.html: 03_clustering.Rmd reports/02_differential_expression_analysis.html
	mkdir -p reports
	Rscript .scripts/compile_rmarkdown.R $< $@

reports/04_downstream_analysis.html: 04_downstream_analysis.Rmd reports/03_clustering.html
	mkdir -p reports
	Rscript .scripts/compile_rmarkdown.R $< $@

reports/%.html: %.Rmd
	mkdir -p .results
	mkdir -p reports
	Rscript .scripts/compile_rmarkdown.R $< $@


results/all_go_terms.tsv: .results/clustering_labels.txt
	Rscript run_pathways.R
	Rscript run_go_terms.R


clean:
	rm *.html
	rm -rf reports
	rm -rf .results
	mkdir -p reports
