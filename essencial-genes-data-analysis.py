'''
++++++ Imports ++++++
'''
import sys
import os
import shutil
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

def getGenesModels(genes, models):
    genes = removeTabs(genes)
    models = removeTabs(models)

    for i in range (len(models)):
        rename = models[i].split('_')
        models[i] = rename[0]

    essencial_models = []
    for i in range (len(models)):
        for j in range (len(genes)):
            if models[i] == genes[j]:
                essencial_models.append(models[i])

    with open('essencial-models-atualizado.txt', 'w') as f:
        f.write('\n'.join(essencial_models))

    df = pd.read_csv(genes_uniparc, sep='\t', header=0)
    df['UniParc'] = df['UniParc'].str.strip()
    filt = df['UniParc'].isin(essencial_models)
    return df[filt]

def stackedBarPlot(data):
    #Set theme
    sns.set_theme(style="whitegrid")
    #plt.rcParams["axes.labelsize"] = 20
    sns.set(font_scale=1.8)

    # Initialize the matplotlib figure
    f, ax = plt.subplots(figsize=(14, 12.75))

    # Load data
    plot = pd.read_excel(data, header=0 ).sort_values('genes', ascending=False)

    # Plot total
    sns.barplot(x='genes', y='pathway', data=plot, label='proteínas', color='#B2DFDB') #linewidth=1.5, edgecolor='.2'
    sns.set_context(rc = {'patch.linewidth': 0.0})

    #sns.set_xlabel("X Label",fontsize=20)
    #sns.set_ylabel("Y Label",fontsize=20)

    # Plot models
    sns.barplot(x='models', y='pathway', data=plot, label='Proteínas com modelo estrutural', color= '#80CBC4')
    sns.set_context(rc = {'patch.linewidth': 0.0})



    # Add a legend and informative axis label
    ax.legend(ncol=2, loc='lower right', frameon=True)
    ax.set(xlim=(0, 12), ylabel='', xlabel='Reações de chokepoint de P. aeruginosa CCBH4851')

    # Remove y axis and save figure
    sns.despine(left=True, bottom=False)
    plt.tight_layout()
    plt.savefig('cp-genes-green3.png', dpi=300)

'''
++++++ File Paths ++++++
'''

genes = '/home/camila/LMDM/P.Aeruginosa/CCBH4851/CCBH4851-metabolic-network/essencial-genes-analysis/essencial-genes-uniparc.txt'
models = '/home/camila/LMDM/P.Aeruginosa/CCBH4851/ThreeFilters-Standalone/good_models.lst'
genes_uniparc = '/home/camila/LMDM/P.Aeruginosa/CCBH4851/CCBH4851-metabolic-network/essencial-genes-analysis/uniprot-essencial-genes-search.csv'
data ='/home/camila/LMDM/P.Aeruginosa/CCBH4851/CCBH4851-metabolic-network/essencial-genes-analysis/bar-plot-traduzido.xlsx'
dataCP= '/home/camila/LMDM/P.Aeruginosa/CCBH4851/CCBH4851-metabolic-network/chokepoint-analysis/CP-analysis/bar-plot-cp-traduzido.xlsx'

file = getGenesModels(genes, models)
file.to_csv('/home/camila/LMDM/P.Aeruginosa/CCBH4851/CCBH4851-metabolic-network/essencial-genes-analysis/essencial-genes-models-atualizado.csv', index = False)
stackedBarPlot(dataCP)
