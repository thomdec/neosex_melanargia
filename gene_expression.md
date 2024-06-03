# Gene expression analysis

## Create the neo-W pseudo-reference 

To minimise reference bias (i.e. poor mapping due to divergence between the query and the reference), we created a pseudo-reference with neo-W alleles (here we used the neo-W alleles from female PT_MI_61, as it was the same sample that was RNAseq'ed).

```bash
bcftools consensus -f melanargia_ines.PT_MI_8.v2_0.fasta -H 1 PT_MI_61.neoW.vcf.gz > PT_MI_61.neoW.fasta
```

## RNAseq alignment

An HISAT2 index must first be build for the reference and the neo-W pseudo-reference:

```bash
hisat2-build -p 10 melanargia_ines.PT_MI_8.v2_0.sequences.red_softmasked.fasta melanargia_ines.PT_MI_8.v2_0.sequences.red_softmasked

hisat2-build -p 10 PT_MI_61.neoW.fasta PT_MI_61.neoW
```

RNAseq data can then be aligned to the reference genome and the neo-W pseudo reference:

```bash
hisat2 -p 10 -x 0_ref/melanargia_ines.PT_MI_8.v2_0.sequences.red_softmasked -1 8_RNAseq/0_trimmed_reads/PT_MI_61.melanargia_ines.concat.1.p.fq.gz -2 8_RNAseq/0_trimmed_reads/PT_MI_61.melanargia_ines.concat.2.p.fq.gz --summary-file 8_RNAseq/1_hisat2/PT_MI_61.v2_0.vs.melanargia_ines.PT_MI_8.v2_0.sequences.red_softmasked.summary.txt | samtools sort -@ 9 -T /scratch/tdecroly/PT_MI_61.v2_0.vs.melanargia_ines.PT_MI_8.v2_0.sequences.red_softmasked.temp.bam -o 8_RNAseq/1_hisat2/PT_MI_61.v2_0.vs.melanargia_ines.PT_MI_8.v2_0.sequences.red_softmasked.bam

hisat2 -p 10 -x 0_ref/PT_MI_61.neoW -1 8_RNAseq/0_trimmed_reads/PT_MI_61.melanargia_ines.concat.1.p.fq.gz -2 8_RNAseq/0_trimmed_reads/PT_MI_61.melanargia_ines.concat.2.p.fq.gz --summary-file 8_RNAseq/1_hisat2/PT_MI_61.vs.PT_MI_61.neoW.txt | samtools sort -@ 9 -T PT_MI_61.vs.PT_MI_61.neoW.temp.bam -o 8_RNAseq/1_hisat2/PT_MI_61.vs.PT_MI_61.neoW.bam
```

## Extract neo-W and neo-Z specific reads

The individual for which RNAseq data is available was also whole-genome sequenced. Neo-W-specific alleles were identified as described in the neosex_phasing section, and encoded in a .tsv file `neoW.target_alleles.tsv`. Same for the neo-Z specific alleles.
Neo-W reads were extracted for the reads aligned to the pseudo-reference

```bash
python extract_reads.py -i 1_hisat2/PT_MI_61.vs.PT_MI_61.neoW.bam -o 2_neoW_bwa_mem/PT_MI_61.vs.PT_MI_61.neoW.neoW_reads.bam -t target_alleles/neoW.target_alleles.tsv

python extract_reads.py -i 8_RNAseq/1_hisat2/PT_MI_61.v2_0.vs.melanargia_ines.PT_MI_8.v2_0.sequences.red_softmasked.bam -o 3_neoZ_bwa_mem/PT_MI_61.v2_0.vs.melanargia_ines.PT_MI_8.v2_0.sequences.red_softmasked.neoZ_reads.bam -t target_alleles/neoZ.target_alleles.tsv
```

The python script is pretty much an unmodified from Simon Martin's `filterSAMbyTargetBase.py` available at https://github.com/simonhmartin/genomics_general/blob/master/SAM_processing/filterSAMbyTargetBase.py I only change one line, to better read in my tsv file with neo-W-specific alleles. 


## Quantifying gene expression with HTseq-count

```bash
samtools sort 2_neoW_bwa_mem/PT_MI_61.vs.PT_MI_61.neoW.neoW_reads.bam | htseq-count -f bam - ../1_clean_vcfs/gene_annotation/melanargia_ines.PT_MI_8.v2_0.sequences.red_repeats.augustus.gt.gff3 -i ID -t gene -r pos > neoW.gene.pos.read_count.txt

samtools sort 3_neoZ_bwa_mem/PT_MI_61.v2_0.vs.melanargia_ines.PT_MI_8.v2_0.sequences.red_softmasked.neoZ_reads.bam | htseq-count -f bam - ../1_clean_vcfs/gene_annotation/melanargia_ines.PT_MI_8.v2_0.sequences.red_repeats.augustus.gt.gff3 -i ID -t gene -r pos > neoZ.gene.pos.read_count.txt
```
