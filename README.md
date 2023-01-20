<p align="center">
  <img src="./img/logo.jpg" width="200">
</p>


`uv` is a small workflow to reliably identify prophages in bacterial genomes. Of course, this is an old problem, but `uv` differs in two regards:

(1) While it uses many different programs and databases, all the bioinformatic nonsense has been taken care of -- parsing, parallelization, more parsing.

(2) `uv` takes a brute force appoach to the prophage finding problem. It identifies putative prophage regions by screening proteins against the [Gut Phage Database](https://www.biorxiv.org/content/10.1101/2020.09.03.280214v1) with about 142,000 metagenome-assembles phages. It then refines candidate regions using [`CheckV`](https://www.biorxiv.org/content/10.1101/2020.05.06.081778v1). Optionally, it reannotates the result using a phage-specific reading frame model.


### Why "uv"?

In the lab, temperate phages (intact prophages) can be induced with UV light, which makes them leave their host's genome and enter in a lytic life cycle. And just as this uncovers (pro)phages using brute physical force, `uv` uncovers them using brute data force.


### Run

```bash
# Get code
git clone https://github.com/phiweger/uv
# Set up environment
conda env create -f uv/env.yml && conda activate uv
# Choose a working directoy
WD=/some/directory
# Get databases and test data -- you need about 24 GB disk space
cp get_db.sh $WD && cd $WD && sh get_db.sh
# Turn on the UV light
nextflow run /path/to/uv/main.nf \
    --standalone \
    --results results \
    --uvdb db \
    --genomes metadata.csv \
    --annotate true \
    --phages phages.fasta \
    -resume
# Add --maxram 8 to limit RAM usage, defaults to 8 GB
# Add -resume (with a single "-"!) to rerun workflow with cached results
# phages.fasta is a collection of phages you'd like to search against.
```


### Prepare a custom phage protein database

This is only necessary if you want to use domain-specific phage proteins to search for these in bacterial/ archaeal genomes. For example, we here use the Gut Phage Database which recruits proteins from ... well, the gut. If you're looking for marine phages, you might want to use known phage proteins from this environment.

```bash
# http://ftp.ebi.ac.uk/pub/databases/metagenomics/genome_sets/gut_phage_database/
wget http://ftp.ebi.ac.uk/pub/databases/metagenomics/genome_sets/gut_phage_database/GPD_proteome.faa
gzip GPD_proteome.faa
mv GPD_proteome.faa.gz 2020-09-16_GPD_proteome.faa.gz

DIR=/path/to/gut_phages
singularity pull --dir $SINGULARITY_CACHEDIR docker://nanozoo/mmseqs2:11.e1a1c--55acb62
singularity exec docker://nanozoo/mmseqs2:11.e1a1c--55acb62 /bin/bash
DB=2020-09-16_GPD_proteome.faa.gz
mmseqs createdb $DB targetDB
mmseqs createindex --threads 40 targetDB tmp
tar zcvf 2020-09-16_GPD_proteome_mmseqs2.tar.gz target*

# create with split index
mmseqs createindex --threads 40 --split 4 targetDB tmp
tar zcvf 2020-09-17_GPD_proteome_mmseqs2_split4.tar.gz target*
```


### More Details

`uv` needs a metadata file with two columns, "name" and "path". "name" is a unique identifier, and "path" is the absolute path to your genome assembly:

```csv
name,path
foo,/path/to/genome_1.fasta
bar,/path/to/genome_2.fasta
```

"name" will appear in the results folder for each bacterium screened. The reason to add this explicit identifier "name" is that humans have just come up with way too many ways to name files, so any parser will break for some people.
