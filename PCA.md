# PCA

```python
import allel
import numpy as np
import pandas as pd

# import vcf 
callset = allel.read_vcf('variant.vcf.gz',  fields=['calldata/GT', 'samples'])

g = allel.GenotypeArray(allel.GenotypeDaskArray(callset['calldata/GT'])) # genotype array
ac = g.count_alleles() # allele count array
samples = callset['samples']

# filter singletons and multi-allelic sites
mask = (ac.max_allele() == 1) & (ac[:, :2].min(axis=1) > 1)
masked_g = g.compress(mask, axis=0)

# create the PCA input (i.e. number of alternative alleles per individual at each site)
gn = masked_g.to_n_alt()

# run the PCA
coords, model = allel.pca(gn, n_components=4, copy=True, scaler='patterson', ploidy=2)

# Coordinates of PC1 and PC2 are recorded in coords[:, 0] and coords[:, 1]. The variance explained by each component is recorded in model.explained_variance_ratio_[n -1], where n is the number of the component.
pca = pd.DataFrame({'ID': samples, 'PC1': coords[:, 0], 'PC2': coords[:, 1], 'PC3': coords[:, 2], 'PC4': coords[:, 3]})
```