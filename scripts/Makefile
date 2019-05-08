ALL_RMD_FILES = $(wildcard *.Rmd)
RMD_FILES = $(filter-out manuscript.Rmd, $(ALL_RMD_FILES))
SUPP_RMD_FILES = $(wildcard miscellaneous/*.Rmd)

HTML_FILES = $(patsubst %.Rmd, reports/%.html, $(RMD_FILES))
SUPP_HTML_FILES = $(patsubst %.Rmd, reports/%.html, $(SUPP_RMD_FILES))


all: html pdf

pdf: manuscript.pdf

html: $(HTML_FILES)

complete: html
	make results/all_go_terms.tsv
	make html

supp: $(SUPP_HTML_FILES)

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
