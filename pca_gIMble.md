# Demographic inference

## PCA

```python
import allel
import numpy as np

# import metadata
metadata = pd.read_csv("../../z_metadata/id_sex_region.csv")

# import vcf 
callset = allel.read_vcf('../../1_clean_vcfs/autosomal/melanargia_ines_preprocess_autosomal_lw.vcf.gz',  fields=['calldata/GT', 'samples'])

g = allel.GenotypeArray(allel.GenotypeDaskArray(callset['calldata/GT'])) # genotype array
ac = g.count_alleles() # allele count array
samples = callset['samples']

# filter singletons and multi-allelic sites
mask = (ac.max_allele() == 1) & (ac[:, :2].min(axis=1) > 1)
mask_g = g.compress(mask, axis=0)

# create the PCA input (i.e. number of alternative alleles per individual at each site)
gn = mask_g.to_n_alt()
gn

# run the PCA
coords, model = allel.pca(gn, n_components=4, copy=True, scaler='patterson', ploidy=2)

# Coordinates of PC1 and PC2 are recorded in coords[:, 0] and coords[:, 1]. The variance explained by each component is recorded in model.explained_variance_ratio_[n -1], where n is the number of the component.
```

## gIMble global demography 

Note: I used gIMble v0.7 (i.e. not the latest version). The commands of the current version differ slightly.

gIMble measure is run to create the gIMble zarr repository. It requires the partition of the genome, encoded in a bed file, the samples and their assigned population (encoded in a `.csv` file), the vcf filtered by gIMble preprocess and the chromosomes analysed, with their respective sizes (encoded in a `.tsv` file).

```bash
gIMble measure -g 02_preprocess/melanargia_ines_15_inds.autosomal.genomefile -b 02_preprocess/autosomes_intergenic.bed.gz -s metadata/melanargia_15.regions.samples.csv -v 02_preprocess/melanargia_ines_15_inds.vcf.gz -z autosomes_intergenic
```

gIMble blocks create block of defined length `-l`

```bash
gIMble blocks -z autosomes_intergenic.z -l 64
```

gIMble info returns basic popgen parameters on the population in the dataset 

```bash
gIMble info -z autosomes_intergenic.z
```

gIMble tally tallies mutations in blocks

```bash
gIMble tally -z autosomes_intergenic.z -d 'blocks' -l tally_blocks_k_2
```

Finally, to run gIMble optimize (to model a global demography), one must first create a configuration file. The configuration file is created with the following command. One must the edit the file manually (e.g. using the `nano` command or `vi` for hardcore users) and add the required info: mutation rate etc. 

```bash
gIMble makeconfig -l <label> -t optimize -m <model>
```

gIMble optimize is finally run as follow:

```bash
gIMble optimize -z autosomes_intergenic.z -c gimble.optimize.DIV.div_model.config.ini -t tally_blocks_k_2
```

Information on the gIMble optimize run can be obtained using gIMble query 

```bash
gIMble query -l optimize/tally_blocks_k_2/div_model -z autosomes_intergenic.z

# One can even generate the demes .yml file with gIMble query

gIMble query -l optimize/tally_blocks_k_2/div_model -z autosomes_intergenic.z --demes
```


