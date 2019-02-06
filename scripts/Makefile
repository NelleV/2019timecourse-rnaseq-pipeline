

PIPELINE_STEPS = step0_preprocessing \
		 step2_normalization \
		 step3_DE_analysis \
		 step4_clustering

all:
	for i in $(PIPELINE_STEPS); do $(MAKE) -C $$i all; done
