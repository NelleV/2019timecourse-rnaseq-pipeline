RSCRIPT ?= Rscript


all: doc check build install test

build:
	$(RSCRIPT) .scripts/build_package.R

install:
	$(RSCRIPT) .scripts/install_package.R


check:
	$(RSCRIPT) .scripts/check_package.R

doc:
	$(RSCRIPT) .scripts/build_documentation.R


test:
	$(RSCRIPT) .scripts/test_package.R
