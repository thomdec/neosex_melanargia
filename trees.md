# Phylogenetic rees

## Create a consensus sequences

Here is a bash script to create consensus fasta sequence for each individual in the VCF. The first step is to create a negative mask encoding missing regions of the genome for the individual (one must not encode missing regions as if they had the reference allele). This is done by subtracting called regions of the individual (from the gIMble multi-callable bed output) from the total  genome bed file (chromosome_1   0   51431728/n chromosome_2  0   40278090 ...). The consensus sequence is then obtained using `bcftools consensus`.

```bash
#!/bin/bash

#create_consensus.sh

ind=$1

echo "Processing $ind"

# create negative mask
grep "$ind" callable.bed | bedtools merge | bedtools subtract -a complete_genome.bed -b - > beds/$ind.negative.bed

# extract consensus sequence
bcftools consensus --sample $ind -m beds/$ind.negative.bed -f ref.fasta variant.vcf.gz -p "$ind " > fasta/$ind.fa
```

The script is iterated through individuals with `while read ind; do bash create_consensus.sh $ind; done < individuals.txt`

### Consensus sequence alignment

Consensus fasta file are merge with `cat fasta/* > merged.fa`.

```bash
# alignment with mafft
mafft --auto merged.fa > merged.mafft.fa

# triming alignments with trimal
trimal -in merged.mafft.fa -out merged.mafft.trimal.fa -gappyout
```

### Maximum likelihood tree construction with IQTREE-2

```bash
iqtree2 -s merged.mafft.trimal.fa -m TEST -B 1000 -T AUTO --prefix merged
```
