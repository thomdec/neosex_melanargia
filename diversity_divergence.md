# Estimation of nucleotide diversity and divergence

## Global estimates of $\pi$ and d_{xy}

Estimation of diversity and divergence are based on allelic difference between samples, as well as the total number of pairwise comparison. 

- The number of allelic difference is computed from the VCF using scikit-allel. 

- The number of pairwise comparison conducted is computed from the multi-callable bed array &mdash; a file containing all the callable position for each individual.

```
# example of a multi-callable bed array
#The array is composed of the chromosome in the first column, the intervals (start, end) in the second and third two columns, as well as a $n \times k$ matrix, where n is the number of intervals, and k is the number of samples. Callable intervals are recorded with $1$, missing intervals are recorded with $0$.
# the order of samples is in the callable matrix is the same as in the VCF file, for ease of computation
contig1 0   100 1   0   1
contig1 100 150 1   1   1
```

Filtered VCF were created with [gIMble repo](https://github.com/DRL/gimbleprep) (see [Variant callaning](variant_calling_and_filtering.md)). One can use the output of gimbleprep to make the multicallable bed array. When gimbleprep is run with the `-k` parameter (to keep the intermediate files), the hidden directory it produces contains:
    - the callable regions for each individual (i.e. regions within the minimum depth set and 3 times the average depth, obtained with mosdepth quantised)
    - a VCF file of failed variants (filtered out in the clean VCF, e.g. low quality sites, indels...)

The multicallable bed array can be created as follow:

```bash
# create the multicallable array
bedtools multiinter -i DZ_MI_1624.callable.bed DZ_MI_1681.callable.bed DZ_MI_1685.callable.bed MA_MI_1620.callable.bed TN_MI_1619.callable.bed ES_MI_1680.callable.bed ES_MI_1682.callable.bed ES_MI_1683.callable.bed ES_MI_1684.callable.bed ES_MI_1686.callable.bed PT_MI_61.callable.bed PT_MI_86.callable.bed PT_MI_7.callable.bed PT_MI_8.callable.bed ES_MI_1647.callable.bed -names DZ_MI_1624 DZ_MI_1681 DZ_MI_1685 MA_MI_1620 TN_MI_1619 ES_MI_1680 ES_MI_1682 ES_MI_1683 ES_MI_1684 ES_MI_1686 PT_MI_61 PT_MI_86 PT_MI_7 PT_MI_8 ES_MI_1647 | cut -f 1-3,6- | gzip > multiinter.callable.array.bed.gz

# remove the SNPs filtered out of the VCF from the multicallable array
bcftools view vcf.filtered.vcf.gz -i "%FILTER!='PASS'" | bcftools query -f '%CHROM\t%POS0\t%END\t%FILTER\n' | gzip > vcf.fails.bed.gz
bedtools subtract -a multiinter.callable.array.bed.gz -b vcf.fails.bed.gz | gzip > multiinter.callable.array.clean.bed.gz
```

Diversity and divergence can then be computed with the following functions.

```python
import allel
import pybedtools
import numpy as np
import pandas as pd
from scipy.special import comb
from itertools import combinations
import subprocess

def diff_within(g, pop):
    """Returns per site allelic differences within a population"""
    ac_pop = g[:, pop].count_alleles()
    an_pop = np.sum(ac_pop, axis=1)
    n_pairs = an_pop * (an_pop - 1) / 2
    n_same = np.sum(ac_pop * (ac_pop - 1) / 2, axis=1)
    n_diff = np.sum(n_pairs - n_same)
    return(n_diff)

def diff_between(g, pop1, pop2):
    """Returns per site allelic differences between populations"""
    ac_pop1 = g[:, pop1].count_alleles(max_allele=3)
    ac_pop2 = g[:, pop2].count_alleles(max_allele=3)
    an_pop1 = np.sum(ac_pop1, axis=1)
    an_pop2 = np.sum(ac_pop2, axis=1)
    n_pairs = an_pop1 * an_pop2 
    n_same = np.sum(ac_pop1 * ac_pop2, axis=1)
    n_diff = np.sum(n_pairs - n_same)
    return(n_diff)

def comp_within(coverage, start_end, pop):
    """XX"""
    sum_cov_pop = np.sum(coverage[:, pop], axis=1)
    n_comp = np.sum(comb(2 * sum_cov_pop, 2) * (start_end[:, 1] - start_end[:, 0]))
    return(n_comp)

def comp_between(coverage, start_end, pop1, pop2):
    """XX"""
    sum_cov_pop1 = np.sum(coverage[:, pop1], axis=1)
    sum_cov_pop2 = np.sum(coverage[:, pop2], axis=1)
    n_comp = np.sum(4 * sum_cov_pop1 * sum_cov_pop2 * (start_end[:, 1] - start_end[:, 0]))
    return(n_comp)

def calc_pi(g, coverage, start_end, pop):
    """Returns pi within a population"""

    n_diff = diff_within(g, pop)
    n_comp = comp_within(coverage, start_end, pop)
    pi = n_diff / n_comp
    return(pi)

def calc_dxy(g, coverage, start_end, pop1, pop2):
    """"Returns dxy between two populations"""

    n_diff = diff_between(g, pop1, pop2)
    n_comp = comp_between(coverage, start_end, pop1, pop2)
    dxy = n_diff / n_comp
    return(dxy)

```

Populations are defined using a .csv file containing the population for each sample. 

```python
pop_df = pd.read_csv(pop_file, names=["id", "pop"])
pops = {}
for pop in unique_pops:
    ids = pop_df[pop_df['pop'] == pop]['id']
    pops[pop] = np.where(np.isin(vcf_header.samples, ids))[0]


# import vcf
vcf_header = allel.read_vcf_headers(vcf_file)
callset = allel.read_vcf(vcf_file, fields=['calldata/GT', 'variants/POS', 'samples'])
pos = callset['variants/POS']
g = allel.GenotypeArray(allel.GenotypeDaskArray(callset['calldata/GT']))

# import multicallable bed file
start_end = np.loadtxt(multicallable_file, usecols=(1, 2), dtype=int)
coverage = np.loadtxt(multicallable_file, usecols=(list(range(3, len(vcf_header.samples) + 3))), dtype=int)

# compute pi
for pop in pops:
    print(pop, calc_pi(g, coverage, start_end, pops[pop]))

# compute dxy
for pop1, pop2 in combinations(pops, 2):
    print(pop1, pop2, calc_dxy(g, coverage, start_end, pops[pop1], pops[pop2]))

# compute heterozygosity for each individual 
for i in range(len(vcf_header.samples)):
    print(vcf_header.samples[i], calc_pi(g, coverage, start_end, [i]))
```

### Including a mask (eg. 0D and 4D sites)

```python
mask_bed = pybedtools.BedTool("mask_bed_file.bed")
multicallable_bed = pybedtools.BedTool("multicallable_bed_file.bed")

# create mask for the genotype array
variant_bed = pybedtools.BedTool.from_dataframe(pd.DataFrame({'chr' : chrom, 'start' : pos - 1, 'end' : pos}))
mask_variant_bed = variant_bed.intersect(mask_bed, c=True)
mask_variant = np.asanyarray(mask_variant_bed)[:, 3].astype(int) > 0

# create mask for the multicallable bed file
multicallable_bed = pybedtools.BedTool(multicallable_file)
mask_multicallable_bed = multicallable_bed.intersect(mask_bed)
mask_multicallable_array = np.asanyarray(mask_multicallable_bed)[:, 1:].astype(int)

# compute pi
for pop in pops:
    print(pop, calc_pi(g[mask_variant], mask_multicallable_array[:, 2:], mask_multicallable_array[:, :2], pops[pop]))

# compute dxy
for pop1, pop2 in combinations(pops, 2):
    print(pop1, pop2, calc_dxy(g[mask_variant], mask_multicallable_array[:, 2:], mask_multicallable_array[:, :2], pops[pop1], pops[pop2]))
```

## Windowed estimates of $\pi$ and d_{xy}

```python
def create_windows(window_size, chr_len):
    """Create windows"""

    window_starts = [*range(0, int(chr_len), window_size)]
    window_ends = [*range(0 + window_size, int(chr_len) + window_size, window_size)]
    window_ends[-1] = int(chr_len)
    return((window_starts, window_ends))
```

```python
# get the chromosomes lengths
data = subprocess.check_output("bcftools index -s " + vcf_file, shell=True, text=True)
chromosomes = pd.DataFrame([x.split('\t') for x in data[:-1].split('\n')], columns=["id", "len", "vcf_line"])

# create an index for the multicallable bed file
bed_chrm = np.loadtxt(multicallable_file, usecols=0, dtype=str)
skiprows, max_rows = list(), list()
for i, chr in chromosomes.iterrows():
    chr_idx = np.where(bed_chrm == chr.id)
    skiprows.append(np.min(chr_idx))
    max_rows.append(np.max(chr_idx) - np.min(chr_idx) + 1)

chromosomes['skiprows'] = skiprows
chromosomes['max_rows'] = max_rows

pi_file = open('output.pi.tsv', 'w')
dxy_file = open('output.dxy.tsv', 'w')
window_size = 10000

pi_file.write('chromosome' + '\t' + 'start' + '\t' +  'end' + '\t' +  'pop' + '\t' +  'pi' + '\n')
dxy_file.write('chromosome' + '\t' + 'start' + '\t' +  'end' + '\t' +  'pop1' + '\t' +  'pop2' + '\t' +  'dxy' + '\n')

for iteration, chr in chromosomes.iterrows():

    callset = allel.read_vcf(vcf_file, region=chr.id, fields=['calldata/GT', 'variants/POS', 'samples'], tabix=True)
    pos = callset['variants/POS']
    samples = callset['samples']
    g = allel.GenotypeArray(allel.GenotypeDaskArray(callset['calldata/GT']))
    bed_array = np.loadtxt(multicallable_file, usecols=(list(range(1, len(samples) + 3))), skiprows=chr.skiprows, max_rows=chr.max_rows, dtype=int)
    window_starts, window_ends = create_windows(window_size, chr.len)

    for i in range(len(window_starts)):
        
        mask_window = (pos >= window_starts[i]) & (pos <= window_ends[i])
        win_multicallable_array = intersect(bed_array, window_starts[i], window_ends[i])

        for pop in pops:
            win_pi = calc_pi(g[mask_window], win_multicallable_array[:, 2:], win_multicallable_array[:, :2], pops[pop])
            pi_file.write(str(chr.id) + '\t' + str(window_starts[i]) + '\t' +  str(window_ends[i]) + '\t' +  pop + '\t' +  str(win_pi) + '\n')

        for pop1, pop2 in combinations(pops, 2):
            win_dxy  = calc_dxy(g[mask_window], win_multicallable_array[:, 2:], win_multicallable_array[:, :2], pops[pop1], pops[pop2])
            dxy_file.write(str(chr.id) + '\t' + str(window_starts[i]) + '\t' +  str(window_ends[i]) + '\t' +  pop1 + '\t' +  pop2 + '\t' +  str(win_dxy) + '\n')

pi_file.close()
dxy_file.close()
```