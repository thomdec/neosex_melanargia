# Single tree test

## Find variant types from phased data

```{python}
import numpy as np
import pandas as pd
from itertools import permutations, combinations
from scipy.special import binom, comb
from scipy.stats import mode
import allel

def fold(genotypes):
    if np.sum(genotypes) > 3:
        mutations = [0 if x == int(mode(genotypes, keepdims=True)[0]) else 1 for x in genotypes]
    else:
        mutations = [0 if x == 0 else 1 for x in genotypes]
    return mutations

def polarize(genotypes):
    outg = genotypes[0]
    ing = genotypes[1:]
    poll = [0 if x == outg else 1 for x in ing]
    return poll

callset = allel.read_vcf('neosex.haploid.no_missing.4D.vcf.gz', fields=['calldata/GT', 'variants/POS', 'samples', 'variants/ALT', 'variants/REF'])
h = allel.HaplotypeArray(callset['calldata/GT'][:, :, 0])
samples = callset['samples']
neoW = np.asarray(range(7))

n = len(neoW)
idx = neoW
gt = h[:, idx]
variants = np.mod(np.sum(gt, axis=1), n) != 0
folded_genotypes = [fold(i) for i in gt[variants]] 
# all variants were bi-allelic. Multi-allelic sites must be filtered out prior to the analysis
mutations, counts = np.unique(folded_genotypes, axis=0, return_counts=True)
# mutation is an array with all the unique variant types (eg: [1, 0, 0, 0, 1, 0, 1]) and their associated count
```

## Expected variant type spectrum

One can find the expected distribution of variant types under a neutral model of evolution with recombination

```{python}
sfs = [1/i for i in range(1, n)]
folded_sfs = np.add(sfs, np.flip(sfs))[:int(np.floor(n/2))]
sfs_norm = folded_sfs/np.sum(folded_sfs)
expected_freq = sfs_norm[np.sum(mutations, axis=1) - 1] / binom(7, np.sum(mutations, axis=1)).astype(int)
```

## Tree from the variant types

A tree was constructed by hand from the variant types (taking into account the most prevalent non-singleton variant types) 
Afterwards, we can find the mutations which are consistent or inconsistent with the majority tree:

```{python}
true_mutations = np.array([[0, 0, 0, 0, 0, 0, 1], 
          [0, 0, 0, 0, 0, 1, 0],
          [0, 0, 0, 0, 1, 0, 0],
          [0, 0, 0, 1, 0, 0, 0],
          [0, 0, 1, 0, 0, 0, 0],
          [0, 1, 0, 0, 0, 0, 0],
          [1, 0, 0, 0, 0, 0, 0],
          [1, 0, 0, 0, 0, 0, 1],
          [1, 0, 0, 0, 1, 0, 1],
          [0, 1, 1, 0, 0, 1, 0],
          [0, 0, 1, 0, 0, 1, 0]])


true_mutation_idx = [int(np.argwhere(np.all(mutations == x, axis=1))) for x in true_mutations]
false_mutation_idx = np.delete(np.arange(len(mutations)), true_mutation_idx)
```