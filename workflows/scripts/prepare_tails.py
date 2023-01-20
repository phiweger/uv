'''
Tail spike proteins:

- https://www.nature.com/articles/s41598-021-81063-4#data-availability
- https://github.com/dimiboeckaerts/BacteriophageHostPrediction/blob/master/RBP_database.csv

Idea:

Train a phage protein recognizer by embedding proteins w/ ESM (Facebook) or some other model and then from this train a VAE and then calculate if the protein is a phage or an outlier (not phage).


```bash
wget https://raw.githubusercontent.com/dimiboeckaerts/BacteriophageHostPrediction/master/RBP_database.csv
```
'''


import screed


with open('RBP_database.csv', 'r') as file, open('tails.faa', 'w+') as out:
    # skip header
    _ = next(file)

    for line in file:
        x = line.strip().split(',')
        name = x[0]
        protein_seq = x[4]

        out.write(f'>{name}\n{protein_seq}\n')
