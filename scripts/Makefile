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

%.html: %.Rmd
	mkdir -p .results
	Rscript .scripts/compile_rmarkdown.R $<


clean:
	rm *.html
	for i in $(PIPELINE_STEPS); do $(MAKE) -C $$i clean; done
