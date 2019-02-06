DATASETS = $(wildcard results/*/config.yml) \
	   $(wildcard results/*/*/config.yml)

OUTPUT = $(patsubst %/config.yml, %/timecourse_de_analysis.tsv, $(DATASETS)) \
	 $(patsubst %/config.yml, %/development_de_analysis.tsv, $(DATASETS)) 

all: $(OUTPUT)

%/timecourse_de_analysis.tsv: run_timecourse_differential_expression_analysis.R
	Rscript run_timecourse_differential_expression_analysis.R $(@D)


%/development_de_analysis.tsv: run_development_differential_expression_analysis.R
	Rscript run_development_differential_expression_analysis.R $(@D)

