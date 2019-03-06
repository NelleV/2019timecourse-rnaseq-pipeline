DATASETS = $(wildcard results/*/config.yml) \
	   $(wildcard results/*/*/config.yml)

OUTPUT = $(patsubst %/config.yml, %/counts_normed.tsv, $(DATASETS))

all: $(OUTPUT)

%/counts_normed.tsv: normalization_upper_quartile.R
	Rscript normalization_upper_quartile.R $(@D)

