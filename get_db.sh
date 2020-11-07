#! /bin/sh

THREADS=8
PROJECT="vn6ru"
CHECKV="https://portal.nersc.gov/CheckV/checkv-db-v0.6.tar.gz"

echo "Downloading gut phage database proteins and PVOG HMMs ..."
osf -p $PROJECT clone db
mv db/osfstorage/* db && rm -r db/osfstorage
echo "Decompress ..."
for i in 0 1 2 3; do unpigz -p8 db/targetDB.idx.${i}; done
mkdir db/phage_proteins
mv db/targetDB* db/phage_proteins

echo "Downloading CheckV database ..."
wget --show-progress -O tmp ${CHECKV}
mkdir db/checkv
echo "Decompress ..."
tar -C db/checkv --strip-components 1 -xf tmp
rm tmp

echo "Prepare test data"
mv db/data.tar.gz .
mkdir data
tar -C data -xf data.tar.gz
rm data.tar.gz
# This is so ugly, I'm sorry:
awk -v pwd="$PWD" 'BEGIN {FS=","; OFS=","; print "name,path"} NR>1 {print $1, pwd"/data/"$2}' db/metadata.csv > metadata.csv
rm db/metadata.csv

echo "Done."