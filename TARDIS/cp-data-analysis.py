'''
++++++ Imports ++++++
'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
'''
++++++ Functions ++++++
'''

def removeTabs(file):
    #Remove white spaces in a file
    with open(file) as f:
        lines = [line.strip() for line in f]
        return lines

def getGeneUniparc(genes, models):
    genes = removeTabs(genes)
    models = removeTabs(models)

    for i in range (len(models)):
        rename = models[i].split('_')
        models[i] = rename[0]

    cp_models = []
    for i in range (len(models)):
        for j in range (len(genes)):
            if models[i] == genes[j]:
                cp_models.append(models[i])

    with open('cp-models.txt', 'w') as f:
        f.write('\n'.join(cp_models))

    df = pd.read_csv(genes_uniparc, sep='\t', header=0)
    df['UniParc'] = df['UniParc'].str.strip()
    filt = df['UniParc'].isin(cp_models)
    return df[filt]

'''
++++++ File Paths ++++++
'''

genes = '/home/camila/LMDM/P.Aeruginosa/CCBH4851/CCBH4851-metabolic-network/chokepoint-analysis/CP-analysis/cp-uniparc-list.txt'
models = '/home/camila/LMDM/P.Aeruginosa/CCBH4851/ThreeFilters-Standalone/good_models.lst'
genes_uniparc = '/home/camila//LMDM/P.Aeruginosa/CCBH4851/CCBH4851-metabolic-network/chokepoint-analysis/CP-analysis/cp-gene-uniparc.csv'
essencial_genes = '/home/camila/LMDM/P.Aeruginosa/CCBH4851/CCBH4851-metabolic-network/essencial-genes-analysis/essencial-genes-uniparc.txt'

file = getGeneUniparc(genes, models)
file.to_csv('/home/camila/LMDM/P.Aeruginosa/CCBH4851/cp-genes-models-atualizado.csv', index = False)

essencial_genes = removeTabs(essencial_genes)
filt = file['UniParc'].isin(essencial_genes)
candidates = file[filt]
candidates.to_csv('/home/camila/LMDM/P.Aeruginosa/CCBH4851/candidatos_atualizado.csv', index = False)

