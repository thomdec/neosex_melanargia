# Neo-sex chromosome phasing

## Identification of neo-W-specific alleles

```python
import allel
import numpy as np
import pandas as pd

# import the neo-sex chromosome VCF file
callset = allel.read_vcf('../1_clean_vcfs/neoZ/melanargia_ines_preprocess_neoZ_lw.vcf.gz', fields=['calldata/GT', 'variants/POS', 'samples', 'variants/ALT', 'variants/REF'])
g = allel.GenotypeArray(allel.GenotypeDaskArray(callset['calldata/GT']))
pos = callset['variants/POS']
alts = callset['variants/ALT']
refs = callset['variants/REF']
samples = callset['samples']

# define populations
metadata = pd.read_csv("../z_metadata/id_sex_region.csv")

# get populations indices in the sample
males = np.where(np.isin(samples, metadata[(metadata['Class'] == 'Male')]['ID']))[0]
females = np.where(np.isin(samples, metadata[(metadata['Class'] == 'Female')]['ID']))[0]
iberia = np.concatenate((males, females))

# phasing of the neo-sex chromosome
mask_all_het = np.all(g[:, females].is_het(), axis=1) # get site which are heterozygous across all females.
female_ac = g[mask_all_het][:, females].count_alleles()[:, 1:]
male_ac = g[mask_all_het][:, males].count_alleles()[:, 1:]
phased_position_idx, phased_alleles_idx = np.where((female_ac == 7) & (male_ac == 0)) # get positions where females have the same 7 alternative alleles, while males have 0. 
phased_alleles_idx = phased_alleles_idx + 1
pos_phased = pos[mask_all_het][phased_position_idx]
g_phased = g[mask_all_het][phased_position_idx][:, iberia]

## Get position and alleles for the phased genotypes
alleles = np.concatenate((np.stack((refs, refs)).T, alts), axis = 1)[:, 1:]
alleles_phased = alleles[mask_all_het][phased_position_idx]

neoW_alleles = alleles_phased[np.arange(len(alleles_phased)), phased_alleles_idx]
neoW_alleles_df = pd.DataFrame({'chr' : 'melanargia_ines.PT_MI_8.chromosome_14', 'pos' : pos_phased, 'allel' : neoW_alleles})

# export to a tsv file
neoW_alleles_df.to_csv('neoW.target_alleles.tsv', index=False, header=False, sep="\t")
```

## Extraction of neo-W-specific reads 

Reads containing neo-W-specific alleles were extracted from bam file using a python script adapted from `filterSAMbyTargetBase.py` by [Simon Martin](https://github.com/simonhmartin/genomics_general/blob/master/SAM_processing/filterSAMbyTargetBase.py). Online one line of the script was changed, to better read in my tsv file with neo-W-specific alleles. 

```bash
while read ind; do; python extract_neoW_reads.py -i 2_bwa_mem/$ind.vs.melanargia_ines.PT_MI_8.v2_0.sequences.MD.bam -o 2.3_neoW_bwa_mem_final/$ind.vs.melanargia_ines.PT_MI_8.v2_0.sequences.neoW.MD.bam -t target_alleles/neoW.target_alleles.tsv; done < iberian_females.txt
```

## Variant calling and filtering

Bam file containing exclusively neo-W reads were called with freebayes

```bash
freebayes -f 0_ref/melanargia_ines.PT_MI_8.v2_0.sequences.fasta --limit-coverage 500 --use-best-n-alleles 8 --no-population-priors --use-mapping-quality --ploidy 1 --haplotype-length -1 --bam 2.3_neoW_bwa_mem_final/PT_MI_61.vs.melanargia_ines.PT_MI_8.v2_0.sequences.neoW.MD.bam --bam 2.3_neoW_bwa_mem_final/PT_MI_86.vs.melanargia_ines.PT_MI_8.v2_0.sequences.neoW.MD.bam --bam 2.3_neoW_bwa_mem_final/ES_MI_1686.vs.melanargia_ines.PT_MI_8.v2_0.sequences.neoW.MD.bam --bam 2.3_neoW_bwa_mem_final/ES_MI_1684.vs.melanargia_ines.PT_MI_8.v2_0.sequences.neoW.MD.bam --bam 2.3_neoW_bwa_mem_final/ES_MI_1683.vs.melanargia_ines.PT_MI_8.v2_0.sequences.neoW.MD.bam --bam 2.3_neoW_bwa_mem_final/ES_MI_1682.vs.melanargia_ines.PT_MI_8.v2_0.sequences.neoW.MD.bam --bam 2.3_neoW_bwa_mem_final/ES_MI_1680.vs.melanargia_ines.PT_MI_8.v2_0.sequences.neoW.MD.bam | gzip -c - > thomas/01_freebayes/melanargia_ines.vs.melanargia_ines.PT_MI_8.neoW_reads.final.vcf.gz && touch freebayes.melanargia_ines.neoW_reads.final.done
```

Filtering was conducted with [gIMble preprocess](https://github.com/DRL/gIMble/blob/master/cli/preprocess.py)

```bash
gIMble preprocess -f 0_ref/melanargia_ines.PT_MI_8.v2_0.sequences.fasta -v thomas/01_freebayes/melanargia_ines.vs.melanargia_ines.PT_MI_8.neoZ_reads.vcf.gz -b /data/lohse/modelgenomland/analyses/melanargia_ines/2.4_neoZ_bwa_mem -g 2 -q 10 -m 3 -M 3 -t 20 -o thomas/02_preprocess/melanargia_ines.neoZ_haplotypes -k
```
