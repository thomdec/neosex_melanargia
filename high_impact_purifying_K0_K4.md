# High impact variants, $\pi_{0} / \pi_{4}$ and $K_{0} / K_{4}$

## Detection of high impact variants with SnpEff 

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

## Quantifying the number of high impact variant per genes

For some reason, scikit allel cannot handle more than 3 alternative alleles. If a site has > 3 alts, they are simply ignored (this is especially annoying when working with indel data). 

To get around this issue, I created myself the type np.array using the bcftools query function: <br />

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
alleles

# PT_MI_61 (the female for which RNAseq data was available) has is the third sample in the vcf (i.e., its index = 2)

effect_alleles = np.asarray([i.split('|', 1)[0] for i in list(ann)])
effects = np.asarray([i.split('|', 2)[1] for i in list(ann)])
impact = np.asarray([i.split('|', 3)[2] for i in list(ann)])

mask_high_impact = impact == "HIGH" 
# one could also mask only premature stop codons, or frameshifts, using mask = effects == "stop_gained" or mask = effects == "frameshift"
pos_high_impact = pos[mask_high_impact] # find positions of high impact mutations
alts_high_impact = alts[mask_high_impact] # find the allele associated with the high impact mutation

# find the high-impact mutation position and indice in the allele count array
pos_idx, high_impact_idx = np.where(alleles[mask_high_impact] == np.tile(effect_alleles[mask_high_impact], (5, 1)).T) 
ac_high_impact = g[mask_high_impact][:, [2]].count_alleles()

# import genes intervals regions
genic_regions = np.asarray(pd.read_csv('neosex.genic_regions.bed', sep="\t", usecols=(1, 2), names=("start", "end")))

# compute the number of high impact mutation per gene
res = list()
for window in genic_regions:
    mask_gene = (pos_high_impact > window[0]) & (pos_high_impact < window[1])
    n_high_impact = np.sum([np.sum(ac_high_impact[mask_gene][high_impact_idx[mask_gene] == i][:, i]) for i in range(1, np.max(high_impact_idx) + 1)])
    idx_high_impact = np.where(ac_high_impact[mask_gene][high_impact_idx[mask_gene] == 1][:, 1])[0]
    gene_effects = effects[mask_high_impact][mask_gene][idx_high_impact]
    res.append((window[0], window[1], n_high_impact, str(gene_effects)))

df = pd.DataFrame(res, columns=['start', 'end', 'high_impact', 'effect'])
```

## $\pi_{0} / \pi_{4}$

The code for computing diversity on the phased neo-W is slightly different than the one used in [Diversity and divergence](diversity_divergence.md). Since phased neo-W are haploid, the number of comparison between samples is not the same. For a single site in which 7 samples are covered, the number of comparison for haploids is ${{7 \choose 2} = 21}$, while for diploids it is ${7 \times 2 \choose 2} = 91$.

```python

import allel
import numpy as np
import pandas as pd
from scipy.special import comb


# import phased VCF and multicallable.bed for the neo-W at 0D and 4D sites (see neosex_phasing.md for the steps leading to the phased neo-W VCF, see diversity_divergence.md for the steps involved in subsetting the partition)
callset_0 = allel.read_vcf('../1_clean_vcfs/neoW/melanargia_ines.neoW.0_fold_degenerate.lw.vcf.gz', fields=['calldata/GT', 'variants/POS', 'samples'])
samples = callset_0['samples']
h_0 = allel.HaplotypeArray(callset_0['calldata/GT'][:, :, 0])
file = '../1_clean_vcfs/neoW/melanargia_ines.neoW.0_fold_degenerate.multi_callable.bed_array.txt.gz'
start_end_0 = np.loadtxt(file, dtype=int, usecols=[0, 1])
coverage_0 = np.loadtxt(file, dtype=int, usecols=list(range(2, len(samples) + 2)))
pos_0 = callset_0['variants/POS']

callset_4 = allel.read_vcf('../1_clean_vcfs/neoW/melanargia_ines.neoW.4_fold_degenerate.lw.vcf.gz', fields=['calldata/GT', 'variants/POS', 'samples'])
h_4 = allel.HaplotypeArray(callset_4['calldata/GT'][:, :, 0])
file = '../1_clean_vcfs/neoW/melanargia_ines.neoW.4_fold_degenerate.multi_callable.bed_array.txt.gz'
start_end_4 = np.loadtxt(file, dtype=int, usecols=[0, 1])
coverage_4 = np.loadtxt(file, dtype=int, usecols=list(range(2, len(samples) + 2)))
pos_4 = callset_4['variants/POS']

def calc_pi(h, coverage, start_end):
    
    ac_pop = h.count_alleles()
    an_pop = np.sum(ac_pop, axis=1)
    n_pairs = an_pop * (an_pop - 1) / 2
    n_same = np.sum(ac_pop * (ac_pop - 1) / 2, axis=1)
    n_diff = n_pairs - n_same
    tot_diff = np.sum(n_diff)
    sum_cov = np.sum(coverage, axis=1)
    comps = np.sum(comb(sum_cov, 2) * (start_end[:, 1] - start_end[:, 0]))
    tot_comps = np.sum(comps)
    pi = tot_diff / tot_comps
    return(pi)

pi0_pi4 = calc_pi(h_0, coverage_0, start_end_0) / calc_pi(h_4, coverage_4, start_end_4)
```



## $K_0 / K_4$

### Merging the phased neo-W VCF and the diploid neosex chromosome VCF 

To polarise neo-W and neo-Z alleles by the African population, I had to merge the haploid neoW, neoZ VCF with the diploid african VCFs. 

I created pseudo-haplods genomes from african samples, as well as Iberian males:

```bash
#!/bin/bash

# create_pseudohaploids.sh

ind=$1

echo "Processing $ind"

bcftools view -H ../1_clean_vcfs/neoZ/melanargia_ines_preprocess_neoZ_lw.vcf.gz | cut -f 1,2,3,4,5,6,7,8,9 | paste - africa_haps/$ind\_1.txt > $ind.tmp_1.txt

bcftools view -H ../1_clean_vcfs/neoZ/melanargia_ines_preprocess_neoZ_lw.vcf.gz | cut -f 1,2,3,4,5,6,7,8,9 | paste - africa_haps/$ind\_2.txt > $ind.tmp_2.txt

bcftools view -s MA_MI_1620 -h ../1_clean_vcfs/neoZ/melanargia_ines_preprocess_neoZ_lw.vcf.gz | cat - $ind.tmp_1.txt | bgzip > africa_vcfs/$ind.phase_1.vcf.gz

bcftools view -s MA_MI_1620 -h ../1_clean_vcfs/neoZ/melanargia_ines_preprocess_neoZ_lw.vcf.gz | cat - $ind.tmp_2.txt | bgzip > africa_vcfs/$ind.phase_2.vcf.gz

bcftools index africa_vcfs/$ind.phase_1.vcf.gz
bcftools index africa_vcfs/$ind.phase_2.vcf.gz

echo -e "$ind\n$ind" > tmp.txt
echo -e "_1\n_2" | paste tmp.txt - -d "" > names.txt

bcftools merge --force-samples africa_vcfs/$ind.phase_1.vcf.gz africa_vcfs/$ind.phase_2.vcf.gz | bcftools reheader -s names.txt | bcftools view -Oz -o africa_vcfs/$ind.haploids.vcf.gz

bcftools index africa_vcfs/$ind.haploids.vcf.gz

rm $ind.tmp_1.txt
rm $ind.tmp_2.txt
```

```bash
while read ind; do bash create_pseudohaploids.sh $ind; done < inds_africa.txt
```

The pseudohaploids genomes are then merged:

```bash
bcftools merge DZ_MI_1624.haploids.vcf.gz DZ_MI_1681.haploids.vcf.gz DZ_MI_1685.haploids.vcf.gz MA_MI_1620.haploids.vcf.gz TN_MI_1619.haploids.vcf.gz -Oz -o neo_sex.african_haploids.vcf.gz
```

The phased neoW and phased neoZ were merged:

```bash
bcftools merge --force-samples -Oz -o neoZ_neoW.lw.vcf.gz melanargia_ines.neoW.final.lw.vcf.gz melanargia_ines.neoZ.lw.vcf.gz

# the samples were renamed with 

bcftools reheader -s new_names.txt neoZ_neoW.lw.vcf.gz melanargia_ines.neoW.final.lw.vcf.gz -Oz -o named_neoZ_neoW.lw.vcf.gz 
```

The haploid VCFs from Iberian males, african samples, neoW and neoZ are then merged:

```bash
bcftools merge named_neoZ_neoW.lw.vcf.gz ../../8_phasing/africa_vcfs/neo_sex.african_haploids.vcf.gz ../../8_phasing/male_vcfs/neo_sex.male_haploids.vcf.gz -Oz -o neosex.haploids.final.vcf.gz
bcftools index neosex.haploids.final.vcf.gz
```

The multi-callable bed file are created using the following commands: 

```bash
while read ind; do zgrep "$ind" melanargia_ines.neoW.final.multi_callable.bed.gz | bedtools merge | gzip > beds/$ind.neoW.bed.gz;zgrep "$ind" melanargia_ines.neoZ.multi_callable.bed.gz | bedtools merge | gzip > beds/$ind.neoZ.bed.gz; done < iberian_females.txt

while read ind; do grep "$ind" ../neoZ/neoZ.bed | bedtools merge | gzip > beds/$ind.bed.gz done < african_and_male.inds.txt

bedtools multiinter -i DZ_MI_1624.bed.gz DZ_MI_1624.bed.gz DZ_MI_1681.bed.gz DZ_MI_1681.bed.gz DZ_MI_1685.bed.gz DZ_MI_1685.bed.gz MA_MI_1620.bed.gz MA_MI_1620.bed.gz TN_MI_1619.bed.gz TN_MI_1619.bed.gz ES_MI_1680.neoW.bed.gz ES_MI_1680.neoZ.bed.gz ES_MI_1682.neoW.bed.gz ES_MI_1682.neoZ.bed.gz ES_MI_1683.neoW.bed.gz ES_MI_1683.neoZ.bed.gz ES_MI_1684.neoW.bed.gz ES_MI_1684.neoZ.bed.gz ES_MI_1686.neoW.bed.gz ES_MI_1686.neoZ.bed.gz PT_MI_61.neoW.bed.gz PT_MI_61.neoZ.bed.gz PT_MI_86.neoW.bed.gz PT_MI_86.neoZ.bed.gz PT_MI_7.bed.gz PT_MI_7.bed.gz PT_MI_8.bed.gz PT_MI_8.bed.gz ES_MI_1647.bed.gz ES_MI_1647.bed.gz -names DZ_MI_1624_1 DZ_MI_1624_2 DZ_MI_1681_1 DZ_MI_1681_2 DZ_MI_1685_1 DZ_MI_1685_2 MA_MI_1620_1 MA_MI_1620_2 TN_MI_1619_1 TN_MI_1619_2 ES_MI_1680_W ES_MI_1680_Z ES_MI_1682_W ES_MI_1682_Z ES_MI_1683_W ES_MI_1683_Z ES_MI_1684_W ES_MI_1684_Z ES_MI_1686_W ES_MI_1686_Z PT_MI_61_W PT_MI_61_Z PT_MI_86_W PT_MI_86_Z  PT_MI_7_1 PT_MI_7_2 PT_MI_8_1 PT_MI_8_2 ES_MI_1647_1 ES_MI_1647_2 | cut -f1-5 | gzip > haploids_neo_sex.multicallable.bed.gz
```

When working with variant-only VCF file, one cannot simply merge VCF file assuming that sites missing in one VCF are indeed missing (they may be called for the reference allele across all individuals...). I thus restricted my analysis to sites covered across all samples so that sites missing in the VCF are not missing but the reference allel.

```bash
zgrep -P "\t30\t" haploids_neo_sex.multicallable.bed.gz | gzip > haploids_neo_sex.conserved_multicallable.bed.gz

bcftools view -R beds/haploids_neo_sex.conserved_multicallable.bed.gz neosex.haploids.final.vcf.gz | bcftools +missing2ref -Oz -o neosex.haploid.no_missing.vcf.gz
```

I follow the procedure described in [Diversity and divergence](diversity_divergence.md) to subset 4D and 0D sites.

### Computing $K_0 / K_4$

```python
import allel
import numpy as np
import pandas as pd

callset_0 = allel.read_vcf('neosex.haploid.no_missing.0D.vcf.gz', fields=['calldata/GT', 'variants/POS', 'samples'])
samples = callset_0['samples']
h_0 = allel.HaplotypeArray(callset_0['calldata/GT'][:, :, 0])
pos_0 = callset_0['variants/POS']
file = 'haploids_neo_sex.0D.conserved_multicallable.bed.gz'
df = pd.read_csv(file, usecols=[1, 2], sep="\t", index_col=False)
start_end_0 = np.asarray(df)
coverage_0 = np.ones((len(start_end_0), len(samples)), dtype=int)

callset_4 = allel.read_vcf('neosex.haploid.no_missing.4D.vcf.gz', fields=['calldata/GT', 'variants/POS', 'samples'])
h_4 = allel.HaplotypeArray(callset_4['calldata/GT'][:, :, 0])
pos_4 = callset_4['variants/POS']
file = 'haploids_neo_sex.4D.conserved_multicallable.bed.gz'
df = pd.read_csv(file, usecols=[1, 2], sep="\t", index_col=False)
start_end_4 = np.asarray(df)
coverage_4 = np.ones((len(start_end_4), len(samples)), dtype=int)

neoW = np.asarray(range(7))
neoZ = np.asarray(range(7, 14))
algeria = np.asarray(range(14, 20))


## Import gene intervals
genic_regions = np.asarray(pd.read_csv('neosex.genic_regions.bed', sep="\t", usecols=(1, 2), names=("start", "end")))


def gene_k(h, pos, coverage, start_end, pop_idx):
    """Returns the number of derived mutations, polarised by the algerian population"""

    mask_covered = h.count_missing(axis=1) == 0
    h_no_missing = h[mask_covered]
    mask_all_callable = np.all(coverage != 0, axis=1)
    ac_algeria = h_no_missing[:, algeria].count_alleles()
    pos_idx, array_idx = np.where(ac_algeria == len(algeria))
    ac_pop = h_no_missing[pos_idx][:, pop_idx].count_alleles(max_allele=3)

    res = list()
    for window in genic_regions:

        window_start = window[0]
        window_end = window[1]

        mask_win = np.any((start_end[mask_all_callable] >= window_start) & (start_end[mask_all_callable] <= window_end), axis=1)
        win_start_end = start_end[mask_all_callable][mask_win]

        if len(win_start_end) == 0:
            k = None

        else:
            mask_window = (pos[mask_covered][pos_idx]  >= window_start) & (pos[mask_covered][pos_idx] <= window_end)
            derived = np.sum((np.sum(ac_pop[mask_window], axis=1) - ac_pop[mask_window][np.arange(len(ac_pop[mask_window])), array_idx[mask_window]]))
            called_sites = np.sum(win_start_end[:, 1] - win_start_end[:, 0])
            k = derived / called_sites


        res.append((window_start, window_end, k))
    
    return(np.asarray(res))

k_0 = gene_k(h_0, pos_0, coverage_0, start_end_0, neoW)
k_4 = gene_k(h_4, pos_4, coverage_4, start_end_4, neoW)

mask = (k_0[:, 2].astype(float) + k_4[:, 2].astype(float)) > 0
df = pd.DataFrame(genic_regions[mask], columns=['start', 'end'])
df["k_0"] = k_0[mask][:, 2].astype(float)
df["k_4"] = k_4[mask][:, 2].astype(float)
df['k0_k4'] = k_0[mask][:, 2].astype(float) / k_4[mask][:, 2].astype(float) 
```