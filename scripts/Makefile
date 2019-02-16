RMD_FILES = $(wildcard *.Rmd)
HTML_FILES = $(patsubst %.Rmd, %.html, $(RMD_FILES))
R_FILES = $(patsubst %.Rmd, tmp/%.R, $(RMD_FILES))

PIPELINE_STEPS = step0_preprocessing \
		 step2_normalization \
		 step3_DE_analysis \
		 step4_clustering

all: html

html: $(HTML_FILES)


pipeline:
	for i in $(PIPELINE_STEPS); do $(MAKE) -C $$i all; done


03_clustering.html: 03_clustering.Rmd 02_differential_expression_analysis.html
	Rscript .scripts/compile_rmarkdown.R $<

04_downstream_analysis.html: 04_downstream_analysis.Rmd 03_clustering.html
	Rscript .scripts/compile_rmarkdown.R $<


%.html: %.Rmd
	mkdir -p .results
	Rscript .scripts/compile_rmarkdown.R $<

clean:
	rm *.html
	for i in $(PIPELINE_STEPS); do $(MAKE) -C $$i clean; done
