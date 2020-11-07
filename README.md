## README

![](img/logo.jpg)



```bash
# Set up environment
conda env create -f env.yml && conda activate uv
# Get code
git clone github.com/phiweger/uv
# Choose a working directoy
$WD=/some/directory
# Get databases and test data
cp get_db.sh $WD && cd $WD && sh get_db.sh
# Turn on the UV light
nextflow run /path/to/uv/main.nf --results results --db db --genomes metadata.csv
```
