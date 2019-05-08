RMD_FILES = $(wildcard *.Rmd)
SUPP_RMD_FILES = $(wildcard miscellaneous/*.Rmd)

HTML_FILES = $(patsubst %.Rmd, reports/%.html, $(RMD_FILES))
SUPP_HTML_FILES = $(patsubst %.Rmd, reports/%.html, $(SUPP_RMD_FILES))


all: html

html: $(HTML_FILES)

complete: html
	make results/all_go_terms.tsv
	rm -f reports/04_downstream_analysis.html
	make html

supp: $(SUPP_HTML_FILES)

reports/02_differential_expression_analysis.html: 02_differential_expression_analysis.Rmd reports/01_quality_control_normalization.html
	mkdir -p reports
	Rscript .scripts/compile_rmarkdown.R $< $@

reports/03_clustering.html: 03_clustering.Rmd reports/02_differential_expression_analysis.html
	mkdir -p reports
	Rscript .scripts/compile_rmarkdown.R $< $@

reports/%.html: %.Rmd
	mkdir -p .results
	mkdir -p $(@D)
	Rscript .scripts/compile_rmarkdown.R $< $@

results/all_go_terms.tsv: .results/clustering_labels.txt
	Rscript run_pathways.R
	Rscript run_go_terms.R


clean:
	rm *.html
	rm -rf reports
	rm -rf .results
	mkdir -p reports
