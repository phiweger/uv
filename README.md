<p align="center">
  <img src="./img/logo.jpg" width="200">
</p>


`uv` is a small workflow to reliably identify prophages in bacterial genomes. Of course, this is an old problem, but `uv` differs in two regards:

(1) While it uses many different programs and databases, all the bioinformatic nonsense has been taken care of -- parsing, parallelization, more parsing.

(2) `uv` takes a brute force appoach to the prophage finding problem. It identifies putative prophage regions by screening proteins against the [Gut Phage Database](https://www.biorxiv.org/content/10.1101/2020.09.03.280214v1) with about 142,000 metagenome-assembles phages. It then refines candidate regions using [`CheckV`](https://www.biorxiv.org/content/10.1101/2020.05.06.081778v1). Optionally, it reannotates the result using a phage-specific reading frame model.


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
nextflow run /path/to/uv/main.nf --results results --uvdb db --genomes metadata.csv --annotate true
# Add --maxram 8 to limit RAM usage, defaults to 8 GB
# Add -resume (with a single "-"!) to rerun workflow with cached results
```


### Why "uv"?

In the lab, temperate phages (intact prophages) can be induced with UV light, which makes them leave their host's genome and enter in a lytic life cycle. And just as this uncovers (pro)phages using brute physical force, `uv` uncovers them using brute data force.


### Details

`uv` needs a metadata file with two columns, "name" and "path". "name" is a unique identifier, and "path" is the absolute path to your genome assembly:

```csv
name,path
foo,/path/to/genome_1.fasta
bar,/path/to/genome_2.fasta
```

"name" will appear in the results folder for each bacterium screened. The reason to add this explicit identifier "name" is that humans have just come up with way too many ways to name files, so any parser will break for some people.