import numpy as np
import pandas as pd
from openTSNE import TSNE
from umap import UMAP


def normalize(v):
    '''
    stackoverflow, 21030391
    '''
    norm = np.linalg.norm(v)
    if norm == 0: 
       return v
    return v / norm



meta = '../outbreak/microreact.csv'
dfmeta = pd.read_csv(meta)



n_components = 2


fp1 = 'signature.csv'
fp2 = 'containment.phages.csv'
df1 = pd.read_csv(fp1)
df2 = pd.read_csv(fp2)

df = pd.merge(df1, df2, on='name')

redux = UMAP(n_components=10)


df.index = df['name']
df.drop(columns='name', inplace=True)


d = {}
for name, v in df.iterrows():
    v_ = [float(i) for i in v]
    vn = normalize(v_)
    d[name] = vn




# redux = UMAP(n_components=n_components)
# projection = redux.fit_transform(list(d.values()))


redux = TSNE(n_components=n_components)
projection = redux.fit(np.array(list(d.values())))


# dfmeta[dfmeta.id == '19142e05-7365-4b55-abcc-9ba0dec235d2'].country__autocolor.item()



with open('embedding.csv', 'w+') as out:
    header = [f'c{str(i+1)}' for i in range(n_components)]
    out.write('name,country,study,' + ','.join(header) + '\n')  # header
    for leaf, v in zip(d.keys(), projection):
        study = dfmeta[dfmeta.id == leaf].study__autocolor.item()
        country = dfmeta[dfmeta.id == leaf].country__autocolor.item()
        out.write(
            f'{leaf},{country},{study},{",".join([str(i) for i in v])}\n')


with open('embedding_raw.csv', 'w+') as out:
    header = [f'c{str(i+1)}' for i in range(n_components)]
    out.write('name,' + ','.join(header) + '\n')  # header
    for leaf, v in zip(d.keys(), d.values()):
        out.write(
            f'{leaf},{",".join([str(i) for i in v])}\n')


'''R
library(readr)
library(ggplot2)

df <- read_csv('embedding.csv')
ggplot(df, aes(x=c1, y=c2, color=country)) + geom_point(size=0.75) + facet_wrap(~study) + theme_classic()
'''


