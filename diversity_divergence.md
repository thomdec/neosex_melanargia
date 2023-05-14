# Estimation of diversity and divergence

## Global estimate of diversity and divergence

Estimation of diversity and divergence relies on counting allelic difference between samples, as well as the total number of pairwise comparison operated. 

- The number of allelic difference is computed from the VCF using scikit-allel. 

- The number of pairwise comparison conducted is computed from the multi-callable bed file &mdash; a file containing for the individual that are called at each sites. 

```
# example of a multi-callable bed file
contig1 0   100 individual1,individual3
contig1 100 150 individual1,individual2,individual3
```

The bed file is first transformed into a numpy array to speed up computations. The numpy array is composed of the intervals (start, end) in the first two columns, as well as a $n \times k$ matrix, where n is the number of intervals, and k is the number of samples. Callable intervals are recorded with $1$, missing intervals are recorded with $0$.

```
# example of a multi-callable bed file
# the order of samples is in the callable matrix is the same as in the VCF file, for ease of computation
0   100 1   0   1
100 150 1   1   1
```

The bed file is transformed to a numpy array with `bed_to_numpy.py` ([bed_to_numpy script](bed_to_numpy.py)). The script integrates multiprocessing (i.e. parralell processing on multiple cores).

```bash
python bed_to_numpy.py -i <multicallable.bed.gz> -o <output.name.txt.gz> -n <number_of_cores> -v <VCF.file.vcf>


# the output will be automatically gzipped if .gz is added at the end of the filename
# the VCF is required to get the sample positions
```

Heterozygosity $\pi$, $d_{xy}$, $d_a$ and $F_{\mathrm{ST}}$ are computed as follow: 

```python
import allel
import numpy as np
import pandas as pd
from scipy.special import comb
from itertools import combinations, product

def calc_pi(pop_idx):
    """Returns pi within a population"""

    # Calculation of the number of allelic difference within a population
    ac_pop = g[:, idx].count_alleles(max_allele=3)
    an_pop = np.sum(ac_pop, axis=1)
    n_pairs = an_pop * (an_pop - 1) / 2
    n_same = np.sum(ac_pop * (ac_pop - 1) / 2, axis=1)
    n_diff = np.sum(n_pairs - n_same)
    # Calculation of the number of pairwise comparisons 
    sum_cov_pop = np.sum(coverage[:, pop_idx], axis=1)
    tot_comps = np.sum(comb(2 * sum_cov_pop, 2) * (start_end[:, 1] - start_end[:, 0]))
    pi = tot_diffs / tot_comps
    return(pi)

def calc_dxy(pop1, pop2):
    """"Returns dxy between two populations"""

    # Computation of the number of allelic difference between two populations
    ac_pop1 = g[:, idx1].count_alleles(max_allele=3)
    ac_pop2 = g[:, idx2].count_alleles(max_allele=3)
    an_pop1 = np.sum(ac_pop1, axis=1)
    an_pop2 = np.sum(ac_pop2, axis=1)
    n_pairs = an_pop1 * an_pop2 
    n_same = np.sum(ac_pop1 * ac_pop2, axis=1)
    n_diff = np.sum(n_pairs - n_same)

    # Computation of the total number of comparisons between the two populations 
    sum_cov_1 = np.sum(coverage[:, pop1], axis=1)
    sum_cov_2 = np.sum(coverage[:, pop2], axis=1)
    tot_comps = np.sum(2 * sum_cov_1 * 2 * sum_cov_2 * (start_end[:, 1] - start_end[:, 0]))
    dxy = tot_diff / tot_comps
    return(dxy)

# the VCF and the multicallable numpy array are imported

callset = allel.read_vcf(<vcf.gz>, fields=['calldata/GT', 'variants/POS', 'samples'])
g = allel.GenotypeArray(allel.GenotypeDaskArray(callset['calldata/GT']))
pos = callset['variants/POS']
samples = callset['samples']

multicallable_file = <multicallable.txt.gz>
start_end = np.loadtxt(multicallable_file, dtype=int, usecols=(0, 1))
coverage = np.loadtxt(multicallable_file, dtype=int, usecols=(list(range(2, len(samples) + 2))))


# To define populations, one can import a .csv, or a .tsv file with the samples names (ID) and their respective populations (Population). 

#The indices of the populations (i.e. where are the samples from which location in the VCF and the multicallable file) are obtained with:

metadata = pd.read_csv(<metadata.csv>)

pop1 = np.where(np.isin(samples, metadata[(metadata['Population'] == 'Pop1')]['ID']))[0]
pop2 = np.where(np.isin(samples, metadata[(metadata['Population'] == 'Pop2')]['ID']))[0]
pop3 = np.where(np.isin(samples, metadata[(metadata['Population'] == 'Pop3')]['ID']))[0]

# a dictionary of population can be created to iterate over parameters function 
pops = {'pop1' : pop1, 'pop2' : pop2, 'pop3' : pop3}


# Finally, one can compute the statistics with:

## pi
pi_pop1 = calc_pi(pop1)
pi_pop2 = calc_pi(pop2)

## d_xy
dxy = calc_dxy(pop1, pop2)

## F_st
mean_pi = (pi_pop1 + pi_pop2) / 2
fst = (dxy - mean_pi) / (dxy + mean_pi) 

## d_a 
### for d_a, pi_ancestral can be chosen (for instance, it can be set to mean diversity, or to diversity from one population, depending on which is the best assumption for ancestral diversity)

da = dxy - mean_pi

## split time estimate using d_a
mutation_rate = <mutation_rate>
split_time = da / (2 * mutation_rate)
```

To automatised computations across populations, one can use the following codes:

### Heterozygosity 

```python
res = list()
for i, sample in enumerate(samples):
    res.append((samples[i], calc_pi([i])))

## One can export the result using pandas
hets = pd.DataFrame(res, columns=['id', 'heterozygosity'])
hets.to_csv('output_file.csv', index=False)
```

### Diversity

```python
res = list()
for pop in pops:
    res.append(pop, calc_pi(pops[pop]))

diversity = pd.DataFrame(res, columns=['pop', 'pi'])
```

### Divergence

```python
res = list()
for pop1, pop2 in combinations(pops, 2):
    dxy = calc_dxy(pops[pop1], pops[pop2])
    pi_1 = calc_pi(pops[pop1])
    pi_2 = calc_pi(pops[pop2])
    mean_pi = (pi_1 + pi_2) / 2
    fst = (dxy - mean_pi) / (dxy + mean_pi)
    da = dxy - mean_pi # this can be modified, depending on which ancestral d_a is chosen
    res.append((pop1, pop2, dxy, fst, da))

divergence = pd.DataFrame(res, columns=['pop1', 'pop2', 'dxy', 'fst', 'da'])
```

### Pairwise divergence 

```python
# divergence between all sample pairs

res = list()
for i, j in combinations(range(len(samples)), 2):
    dxy = calc_dxy([i], [j])
    res.append((samples[i], samples[j], dxy))

pairwise_div = pd.DataFrame(res, columns=['id1', 'id2', 'dxy'])


# divergence between all pairs from two different populations (i.e. no pair within population)

res = list()
for i, j in product(pop1, pop2):
    dxy = calc_dxy([i], [j])
    res.append((samples[i], samples[j], dxy))

pairwise_div = pd.DataFrame(res, columns=['id1', 'id2', 'dxy'])
```

## Windowed estimate of diversity and divergence


```python
def diff_within(idx):
    """Returns per site allelic differences within a population"""
    ac_pop = g[:, idx].count_alleles()
    an_pop = np.sum(ac_pop, axis=1)
    n_pairs = an_pop * (an_pop - 1) / 2
    n_same = np.sum(ac_pop * (ac_pop - 1) / 2, axis=1)
    n_diff = n_pairs - n_same
    return(n_diff)

def diff_between(idx1, idx2):
    """Returns per site allelic differences between populations"""
    ac_pop1 = g[:, idx1].count_alleles()
    ac_pop2 = g[:, idx2].count_alleles()
    an_pop1 = np.sum(ac_pop1, axis=1)
    an_pop2 = np.sum(ac_pop2, axis=1)
    n_pairs = an_pop1 * an_pop2 
    n_same = np.sum(ac_pop1 * ac_pop2, axis=1)
    n_diff = n_pairs - n_same
    return(n_diff)

def windowed_diff(window_size, pop1, pop2):
    """Returns window mid points, diversity in pop1 and pop2, dxy and fst between pop1 and pop2"""
    interval_start, interval_end = 1, 13036179 # chromosome size

    window_pos_1_list = [*range(interval_start, int(interval_end), window_size)]
    window_pos_2_list = [*range(interval_start + (window_size -1), int(interval_end) + window_size, window_size)]
    window_pos_2_list[-1] = interval_end

    window_list = [list(a) for a in zip(window_pos_1_list,window_pos_2_list)]

    diff_pops = diff_between(pop1, pop2)
    diff_pop1 = diff_within(pop1)
    diff_pop2 = diff_within(pop2)
    cov_pop1 = np.sum(coverage[:, pop1], axis=1)
    cov_pop2 = np.sum(coverage[:, pop2], axis=1) 
    callable_starts = start_end[:, 0]
    callable_ends = start_end[:, 1]

    res = list()
    for window in window_list:

        mid = (window[0] + window[1]) / 2
        mask_window = (pos >= window[0]) & (pos <= window[1])

        win_callable_starts = callable_starts[(callable_starts >= window[0]) & (callable_starts <= window[1])]
        win_callable_ends = callable_ends[(callable_ends >= window[0]) & (callable_ends <= window[1])]

        idx_start_win = int(np.where(callable_starts == win_callable_starts[0])[0])
        idx_end_win = int(np.where(callable_ends == win_callable_ends[-1])[0])

        if win_callable_starts[0] >= win_callable_ends[0]:
            win_callable_starts = np.insert(win_callable_starts, 0, window[0])
            idx_start_win = idx_start_win - 1
        if win_callable_starts[-1] >= win_callable_ends[-1]:
            win_callable_ends = np.append(win_callable_ends, window[1])
            idx_end_win = idx_end_win + 1

        ## pi
        win_diff_pop1 = np.sum(diff_pop1[mask_window])
        win_diff_pop2 = np.sum(diff_pop2[mask_window])
        win_cov_pop1 = cov_pop1[idx_start_win:idx_end_win + 1]
        win_cov_pop2 = cov_pop2[idx_start_win:idx_end_win + 1]
        win_pi_pop1 = win_diff_pop1 / (np.sum(comb(2 * win_cov_pop1, 2) * (win_callable_ends - win_callable_starts)))
        win_pi_pop2 = win_diff_pop2 / (np.sum(comb(2 * win_cov_pop2, 2) * (win_callable_ends - win_callable_starts)))
        win_mean_pi = (win_pi_pop1 + win_pi_pop2) / 2

        ## dxy and fst
        win_diff_between = np.sum(diff_pops[mask_window])
        win_comps = np.sum(2 * win_cov_pop1 * 2 * win_cov_pop2 * (win_callable_ends - win_callable_starts))
        win_dxy = win_diff_between / win_comps
        win_fst = (win_dxy - win_mean_pi) / (win_dxy + win_mean_pi)

        res.append((mid, win_dxy, win_fst, win_pi_pop1, win_pi_pop2))
    
    df = pd.DataFrame(res, columns=["mid", "dxy", "fst", "pi_1", "pi_2"])
    return(df)
```

One can compute diversity and divergence in windows with `df = windowed_diff(20_000, pop1, pop2)`. The resulting dataframe can easily be plotted with matplotlib.

```python
import matplotlib.pyplot as plt

fig, ax = plt.subplots(ncols=1, figsize=(16, 4))
plt.plot(df["mid"], df["dxy"])
plt.xlabel("Position")
plt.ylabel(r'$d_{xy}$')
plt.title(r"Neo-sex chromosome windowed $d_{xy}$")
plt.show()
```

## Diversity and divergence per partition

To compute diversity and divergence by partition, one must first create VCFs and multicallable bed files for each partition. 

I created region bed files &mdash; 0D, 4D, intergenic, exonic, intronic sites &mdash; to subset partitions from the all-sites VCF and multi-callable bed file. 

### Intergenic, exonic and intronic sites

I subseted genic, exonic and intronic regions of autosomes from the annotation bed file. 

```python 
#melanargia_ines.PT_MI_8.v2_0.sequences.red_repeats.augustus.gt.bed.gz is the genome annotation bed file
df = pd.read_csv('melanargia_ines.PT_MI_8.v2_0.sequences.red_repeats.augustus.gt.bed.gz', sep='\t', usecols=[0, 1, 2, 7], names=['chr', 'start', 'end', 'type'])

autosomal_exon = df[(df['type'] == 'exon') & (df['chr'].str.contains('melanargia_ines.PT_MI_8.chromosome_[1-9]\Z|melanargia_ines.PT_MI_8.chromosome_1[0-2]\Z'))]

autosomal_exon[["chr", "start", "end"]].to_csv('autosomes.exon.bed.gz', sep="\t", header=False, index=False)

# the same code was used for genic and intronic regions
```

### Site degeneracy

Site degeneracy was inferred using `codingSiteTypes.py` available at https://github.com/simonhmartin/genomics_general. 

```bash
import pandas

python 4-fold/genomics_general/codingSiteTypes.py --ignoreConflicts -a melanargia_ines.PT_MI_8.v2_0.sequences.red_repeats.augustus.gt.gff3 -f gff3 -r melanargia_ines.PT_MI_8.v2_0.sequences.red_softmasked.fasta -o melanargia_ines.sites_types.txt

#melanargia_ines.PT_MI_8.v2_0.sequences.red_repeats.augustus.gt.gff3 is the genome annotation gff3 file
```

I created regions bed file from the `melanargia_ines.sites_types.txt` output in python. 

```python
import pandas

autosomes = site_types[site_types['scaffold'].str.contains('melanargia_ines.PT_MI_8.chromosome_[1-9]\Z|melanargia_ines.PT_MI_8.chromosome_1[0-2]\Z')]

autosomes_4D = autosomes[(autosomes["degeneracy"] == 4)]
autosomes_4D["start"] = autosomes_4D["position"] - 1

# write bed file
autosomes_4D[["scaffold", "start", "position"]].to_csv('autosomes.4D.bed.gz', sep="\t", header=False, index=False)

```

### Subsetting the all-sites VCF and multi-callable bed file

```bash
# subset the VCF
bcftools view -R <region_file.bed> <all_site.vcf.gz> -Oz -o <output.region_subset.vcf.gz>

# subset multi-callable bed
bedtools intersect -a <all_sites.multicallable.bed.gz> -b <region_file.bed> | bedtools sort | gzip > <output.region_subset.multicallabe.bed.gz>

# the multi-callable bed is then transformed to the numpy array 
python bed_to_numpy.py -i <output.region_subset.multicallabe.bed.gz> -o <output.region_subset.name.txt.gz> -n <number_of_cores> -v <output.region_subset.vcf.gz>
```

To compute diversity and divergence on partitions, one can use the function described above, provided that the appropriate VCF and multi-callable bed file are imported. 


