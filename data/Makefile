

all: shoemaker2015/rawdata

shoemaker2015/rawdata:
	bash sync_shoemaker2015.sh
	mkdir -p shoemaker2015/rawdata
	mv GSE63786_RAW.tar shoemaker2015/rawdata/
	cd shoemaker2015/rawdata/ && tar xvf GSE63786_RAW.tar
