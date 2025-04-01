"""
This program identifies hypermutable coding sequences, which contain higher than expected SNP densities on the
assumption that SNPs are randomly distributed across the entire genome. It relies on the 'Export SNPs' file produced via
Mauve along with GenBank files as reference for each sequence. Use the same GenBank files for both aligning sequences
and as reference.
"""

import numpy as np
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt
from scipy.stats import poisson


# collect and format data
def get_data(bin, file):
    with (open(file, 'r')) as f:
        file = np.loadtxt(f, delimiter='\t', dtype=str)
        file = np.delete(file, 0, 0)
        columns = file.shape[1]
        rm_list = [n for n in range(columns) if n % 3 != 0]
        data = pd.DataFrame(file).drop(rm_list, axis=1)
        snps = np.array(data[0], dtype=str)
        data = pd.DataFrame(data.drop(0, axis=1))
    snp_filter = np.zeros(shape=data.shape, dtype=int)
    for snp, x in zip(snps, range(data.shape[0])):
        count = Counter(snp).most_common(1)[0][0]
        for base, y in zip(snp, range(data.shape[1])):
            if base == count:
                snp_filter[x, y] += 0
            else:
                snp_filter[x, y] += 1
    data = data.mul(snp_filter)
    data.replace('', 0, inplace=True)
    bin_size = bin
    data = data.astype(int).div(bin_size)
    data = np.ceil(data).astype(int)
    maximum = max(list(data.max()))
    column_names = list(data.columns)
    density = np.zeros(shape=(data.shape[1], maximum), dtype=int)
    for x, i in zip(column_names, range(6)):
        for y in data[x]:
            if y != 0:
                density[i, y-1] += 1
    density = density[1:]
    return density, data, maximum


# generate SNP density heatmap
def heatmap():
    window = np.ones(window_size) / window_size
    density_rolling = np.array([np.convolve(d, window, mode='same') for d in density])
    plt.figure(figsize=(12, 8))
    labels = ['SPARK_516_C1', 'SPARK_1954_C1', 'SPARK_1059_C1', 'SPARK_2540_C2', 'SPARK_1507_C1', 'SPARK_1852_C1']
    hm = plt.imshow(density_rolling, cmap='magma', interpolation='nearest', aspect='auto')
    plt.yticks(ticks=range(len(labels)), labels=labels, fontsize=9)
    plt.title('Pairwise Mauve SNP Density', fontsize=18)
    cb = plt.colorbar(hm)
    cb.set_label('Density / ' + str(bin_size) + ' bases', labelpad=10, fontsize=11)
    plt.xlabel('Genome Location (' + str(bin_size) + ') bases', fontsize=11)
    plt.show()


# poisson PMF to determine significance of SNP distribution
def plot_freq(n):
    arr = density[n]
    mean = np.mean(arr)
    max_value = np.max(arr)
    frequency = np.zeros(shape=max_value)
    for x in arr:
        if x != 0:
            frequency[x - 1] += 1
    rolling_avg = np.convolve(frequency, np.ones(window_size) / window_size)
    k = np.arange(max_value)
    pmf = poisson.pmf(k, mean)
    pmf = pmf * len(arr)
    title_dict = {0: 'SPARK_516_C1', 1: 'SPARK_1954_C1', 2: 'SPARK_1059_C1',3: 'SPARK_2540_C2', 4: 'SPARK_1507_C1', 5: 'SPARK_1852_C1'}
    fig, ax = plt.subplots(figsize=(8.5, 11/3))
    ax.bar(np.arange(len(frequency)), frequency, label="Frequency", color='lightblue')
    ax.plot(np.arange(window_size - 1, len(frequency)) - int(window_size / 2), rolling_avg, label="Rolling Average", color='orange')
    ax.plot(k, pmf, label="PMF (Poisson)", color='green', linestyle='--')
    ax.set_xlabel("SNPs / Block")
    ax.set_ylabel("Frequency")
    ax.set_title(f"SNP Density Distribution ({title_dict[n]})")
    ax.legend(loc='best')
    ax.grid(True)
    plt.show()


# # # Global Variables # # #
bin_size = 10000
window_size = 10
snp_files = ['mauve_snps']
arrays = get_data(bin_size, snp_files[0])
# # # Objects for Functions # # #
density = arrays[0]
data = arrays[1]
maximum = arrays[2]
# # # Function Call # # #
heatmap()
# heatmap() will output a heatmap of SNP density - bin_size = 10000; window_size = 10
# plot_freq(n) will output SNP density frequency distributions mapped against the poisson distribution
