"""
Make gene presence/absence matrix (figure 6)
Hypothetical proteins should be BLASTp-ed, this will identify the majority of the hypotheticals
Gene names should be manually inspected and cleaned to remove redundancy
Requires grid overlay annotation
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import ListedColormap

with open('top_20_1b.txt', 'r') as f:
    lines = f.read().splitlines()
    data = []
    for x in lines:
        gene = x.split('\t')
        genome = gene[0]
        gene_name = gene[1]
        if gene_name != "hypothetical protein":
            data.append((genome, gene_name))

df = pd.DataFrame(data, columns=["Genome", "Gene"])
presence_absence = pd.crosstab(df["Gene"], df["Genome"])
labels = ['SPARK_516_C1', 'SPARK_1954_C1', 'SPARK_1059_C1', 'SPARK_2540_C2', 'SPARK_1507_C1', 'SPARK_1852_C1']
cmap = ListedColormap(['tomato', (0.49, 0.99, 0.5)])
plt.figure(figsize=(12, 8))
plt.imshow(presence_absence, cmap=cmap, aspect="auto", interpolation='nearest')
plt.xticks(ticks=np.arange(6), labels=labels)
plt.yticks(ticks=np.arange(len(presence_absence.index)), labels=presence_absence.index, fontsize=6)
plt.title("Gene Presence-Absence Matrix")
plt.show()
