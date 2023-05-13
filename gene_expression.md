# Gene expression analysis

## RNAseq alignment

An HISAT2 index must first be build:

```bash
hisat2-build -p 10 melanargia_ines.PT_MI_8.v2_0.sequences.red_softmasked.fasta melanargia_ines.PT_MI_8.v2_0.sequences.red_softmasked
# melanargia_ines.PT_MI_8.v2_0.sequences.red_softmasked.fasta is the softmasked reference fasta
```

RNAseq data can then be aligned to the reference genome:

```bash
hisat2 -p 10 -x 0_ref/melanargia_ines.PT_MI_8.v2_0.sequences.red_softmasked -1 8_RNAseq/0_trimmed_reads/PT_MI_61.melanargia_ines.concat.1.p.fq.gz -2 8_RNAseq/0_trimmed_reads/PT_MI_61.melanargia_ines.concat.2.p.fq.gz --summary-file 8_RNAseq/1_hisat2/PT_MI_61.v2_0.vs.melanargia_ines.PT_MI_8.v2_0.sequences.red_softmasked.summary.txt | samtools sort -@ 9 -T /scratch/tdecroly/PT_MI_61.v2_0.vs.melanargia_ines.PT_MI_8.v2_0.sequences.red_softmasked.temp.bam -o 8_RNAseq/1_hisat2/PT_MI_61.v2_0.vs.melanargia_ines.PT_MI_8.v2_0.sequences.red_softmasked.bam
```

## Extract neo-W specific reads

The individual for which RNAseq data is available was also whole-genome sequenced. Neo-W-specific alleles were thus identified as described in the neosex_phasing section, and encoded in a .tsv file `neoW.target_alleles.tsv`.

```bash
python extract_neoW_reads.py -i 1_hisat2/PT_MI_61.v2_0.vs.melanargia_ines.PT_MI_8.v2_0.sequences.red_softmasked.bam -o 2_neoW_bwa_mem/PT_MI_61.vs.melanargia_ines.PT_MI_8.v2_0.sequences.red_softmasked.neoW.bam -t target_alleles/neoW.target_alleles.tsv
```

The python script is pretty much unmodified from Simon Martin's `filterSAMbyTargetBase.py` available at https://github.com/simonhmartin/genomics_general/blob/master/SAM_processing/filterSAMbyTargetBase.py I only change one line, to better read in my tsv file with neo-W-specific alleles. 


## Quantifying gene expression with HTseq-count

```bash
samtools sort 2_neoW_bwa_mem/PT_MI_61.vs.melanargia_ines.PT_MI_8.v2_0.sequences.red_softmasked.neoW.bam | htseq-count -f bam - ../1_clean_vcfs/gene_annotation/melanargia_ines.PT_MI_8.v2_0.sequences.red_repeats.augustus.gt.gff3 -i ID -t gene -r pos > neoW.gene.pos.read_count.txt
```
