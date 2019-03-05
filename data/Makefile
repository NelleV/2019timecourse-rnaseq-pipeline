

all: data/shoemaker2015/rawdata

data/shoemaker2015/rawdata:
	bash sync_shoemaker2015.sh
	mkdir -p data/shoemaker2015/rawdata
	mv GSE63786_RAW.tar data/shoemaker2015/rawdata/
	cd data/shoemaker2015/rawdata/ && tar xvf GSE63786_RAW.tar
