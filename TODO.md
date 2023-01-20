## TODO

- [ ] pharokka
- [ ] https://phrogs.lmge.uca.fr/

https://github.com/gbouras13/pharokka#database-installation

- [ ] tail spike database
- [ ] sanitize genomes module
- [ ] fold and search fold database (of phage spike proteins?)
- [ ] reverse containment (see journal, 2022-03-30, `reverse_gather.py`)
- [ ] circular binary segmentation (see methylation work and py code)

---

maybe use shapemers? (cluster phage fee or whatever to discover new ones)
 
https://pubmed.ncbi.nlm.nih.gov/33381814/
https://www.biorxiv.org/content/10.1101/2022.10.11.511548v1
https://www.nature.com/articles/s41594-022-00849-w

could implement this in faltwerk

https://github.com/TurtleTools/afdb-shapemer-darkness/blob/main/scripts/make_shapemers.py

or ask q like: increased resistance means more/ less phages? are they lost if res
increased or do they bring res?

https://www.nature.com/articles/s41564-022-01263-0

not one word of phages!

---

expose mmseqs params in annotation module (search carful, search feet)

---

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC135240/

phage proteomic tree -- classify the (pro)phage!

https://ggdc.dsmz.de/victor.php
https://academic.oup.com/bioinformatics/article/33/21/3396/3933260

---

extend data sources

phage feet (see journal)

https://mjoh223.github.io/jbd-lab.github.io/static/pdf/publications/soto-perez_2019.pdf
https://github.com/jbisanz/HuVirDB/blob/master/readme.md


---

- don't hard code params mmseqs or at least put lower vals
- another segmentation: circular binary ... -- benchmark!
- annotate pVOGs further https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8442406/

http://ftp.ebi.ac.uk/pub/databases/metagenomics/genome_sets/gut_phage_database/README.txt

- PCs_GPD.txt: GPD proteome clustered. Each line represents a protein cluster.
- GPD_proteome_orthology_assignment.txt.gz: Compressed file containing the functional annotation of GPD proteome.

- concat the queries to mmseqs so we don't load the db a million times
- search crispr takes too much memory with metagraph -- any solution? (split index?)
- predict small ORFs 
    - https://www.sciencedirect.com/science/article/pii/S1931312820306193
    - https://academic.oup.com/nar/article/48/3/1029/5556081
- add other annotation resources:
    - phANNs .. https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007845
    - https://github.com/Adrian-Cantu/PhANNs/tree/master/model_training
    - multiphate? .. https://academic.oup.com/bioinformatics/article/35/21/4402/5488969?login=true

- vis raw signal after protein search


```
touch annotation.gff
find GPD_annotations -name '*.gff' | awk '/##FASTA/ {exit} {print}' >> annotation.gff

But code does not match GPD proteome code ...
```


We can even include non-intact prophages bc/ at some point they got in, nevermind their fate in the genome, though of course the change of an intact receptor is smaller.



Search strategy:

- Manually collect receptors, sequences and 3D
- search for more (snowball)
- Where are they? (end of the protein? search beginning and annotate as putative or fold and see if they have some characteristics)

3D, search upper part homology, fold lower part and search, get upper part homology, ...
