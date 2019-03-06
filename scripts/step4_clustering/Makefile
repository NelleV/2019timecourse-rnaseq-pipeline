DATASETS = $(wildcard results/*/config.yml) \
	   $(wildcard results/*/*/config.yml)

OUTPUT = $(patsubst %/config.yml, %/timecourse_de_analysis.tsv, $(DATASETS)) \

all: $(OUTPUT)

%/.tsv: clustering_kmeans.R
	Rscript clustering_kmeans.R $(@D)

