# Loss-of-function variants

## Detection of "high impact" variants with SnpEff 

To use SnpEff with non-model organisms, one must first build a database for the species. 

This can be done by create a directory with the genome name in the SnpEff directory. Inside the directory, there must be the assembly fasta files named `sequences.fa` and the gff3 file named `genes.gff`. The cds fasta must also be present as `cds.fa`, as well as the protein fasta as `protein.fa` (these are to check that the database was built correctly, but they are not stricly required. Checks can be removed using the right argument in the command line). Lastly, the genome name and directory must be added to the snpEff.config file 

```bash
name=melanargia_ines.PT_MI_8.v2_0.sequences.red_softmasked
mkdir $name
echo "$name.genome : $name" >> snpEff.config
```

The database can then be build with:

```bash
java -jar snpEff.jar build -gff3 -v <directory_name> #name of the data directory
```

Genome annotation is run with:

```bash
java -Xmx4g -jar snpEff.jar ann -c snpEff.config -v melanargia_ines.PT_MI_8.v2_0.sequences.red_softmasked -chr melanargia_ines.PT_MI_8.chromosome_14 melanargia_ines.neoZ.snps_structural.vcf.gz | gzip > annotations/melanargia_ines_preprocess.neoZ.snps_structural.ann.vcf.gz	
```

## Detect fixed loss-of-function mutations

For some reason, scikit allel cannot handle more than 3 alternative alleles. If a site has > 3 alts, they are simply ignored (this is especially annoying when working with indel data). 

To get around this issue, one can create the "type" np.array using the bcftools query function: <br />

`bcftools query -f '%POS\t%ALT\n' melanargia_ines.neoW.final.snps_structural.ann.lw.vcf.gz | gzip > melanargia_ines_preprocess.neoW.final.snp_indels.alts.txt.gz`

```python

import allel
import numpy as np
import pandas as pd
from scipy.special import comb

# import vcf and bed files
callset = allel.read_vcf('melanargia_ines.neoW.final.snps_structural.ann.lw.vcf.gz', fields=['calldata/GT', 'variants/POS', 'samples', 'variants/ANN', 'variants/REF', 'variants/ALT'])
g = allel.GenotypeArray(allel.GenotypeDaskArray(callset['calldata/GT']))
pos = callset['variants/POS']
samples = callset['samples']
ann = callset['variants/ANN']
ann = ann.astype(str)
alts = callset['variants/ALT']
refs = callset['variants/REF']


types_df = pd.read_csv('melanargia_ines_preprocess.neoW.final.snp_indels.alts.txt.gz', sep='\t', usecols=[1], names=['types'])
alts = np.full((len(types_df['types']), 4), '', dtype='<U50')

for row, type in enumerate(list(types_df['types'])):
    site_type = type.split(',')
    for i in range(4):
        if len(site_type) > i:
            alts[row, i] = site_type[i]

alleles = np.concatenate((np.tile(refs, (2, 1)).T, alts), axis = 1)[:, 1:]

# PT_MI_61 (the female for which RNAseq data was available) has is the third sample in the vcf (i.e., its index = 2)

effect_alleles = np.asarray([i.split('|', 1)[0] for i in list(ann)])
effects = np.asarray([i.split('|', 2)[1] for i in list(ann)])
impact = np.asarray([i.split('|', 3)[2] for i in list(ann)])

# create high impact variants allele counts
high_impacts = impact == "HIGH"
ac = g[:, iberian_females][high_impacts].count_alleles()

# find the index of high impact alleles
high_impact_idx = np.where(alleles[high_impacts] == np.tile(effect_alleles[high_impacts], (5, 1)).T)[1]

# find fixed high impact alleles in the population
pos_idx, fixed_variants = np.where(ac == 10)
fixed_high_impact = np.equal(fixed_variants, high_impact_idx[pos_idx]) # the index of high impact allel equals the index of the fixed allel

# filter nonsense mutations: premature stop codons, frameshift variants and start loss.
masked_effects = effects[high_impacts][pos_idx][fixed_high_impact]
mask_lof = (masked_effects == "frameshift_variant") | (masked_effects == "start_lost") | (masked_effects == "stop_gained")

pos_lof = pos[high_impacts][pos_idx][fixed_high_impact][mask_lof]
```

## Quantify number of genes with fixed loss-of-function variants

Intersect the position of fixed LOF variants with the position of each genes

```{python}
# a bed file of genes can easily be created from the annotation.gff3
genic_regions = np.asarray(pd.read_csv('genes.bed', sep="\t", usecols=(1, 2), names=("start", "end")))

res = list()
for gene in genic_regions:
    res.append((gene[0], gene[1], np.sum((pos_nonsense > gene[0]) & (pos_nonsense < gene[1]))))

lof_gene_df = pd.DataFrame(res, columns=["start", "end", "n_nonsense"])
```

This process was repeated for: Iberian females, Iberian males and North-African samples (when indexing the genotype array to make the allele count array with .count_alleles())