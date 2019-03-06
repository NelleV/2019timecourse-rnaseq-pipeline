DATASETS = results/varoquaux2019/leaf/meta.tsv \
	   results/varoquaux2019/root/meta.tsv


all: $(DATASETS)


results/varoquaux2019/%/meta.tsv: process_epicon_%.R data/varoquaux2019/%/meta.tsv
	mkdir -p $(@D)
	Rscript process_epicon_$*.R

clean:
	rm -rf results
