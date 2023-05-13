# Phylogenetic trees

## Mitochondrial tree

### Created consensus sequences

Here is a bash script to create consensus fasta sequence for each individual in the mitochondrial VCF. The first step is to create a negative mask encoding missing regions of the genome for the individual. This is done because one must not encode missing regions as if they had the reference allele. This is done by subtracting called regions of the individual (from the gIMble multi-callable bed output) from the total mitochondrial genome bed file (contig1   0   15371). The consensus sequence is then obtained using `bcftools consensus`.

```bash
#!/bin/bash

#create_consensus.sh

ind=$1

echo "Processing $ind"

# create negative mask
grep "$ind" melanargia_ines_mitochondrial.bed | bedtools merge | bedtools subtract -a complete_mitochondria.bed -b - > beds/$ind.negative.bed

# extract consensus sequence
bcftools consensus --sample $ind -m beds/$ind.negative.bed -f Circularized_assembly_1_MI_8.fasta melanargia_ines_mitochondrial.vcf.gz -p "$ind " > fasta/$ind.fa
```

The script is iterated through individuals with `while read ind; do bash create_consensus.sh $ind; done < individuals.txt`


### Consensus sequence alignment

Consensus fasta file are merge with `cat fastas/* > mitochondria.iberian_females.fa`, the `fastas` directory containing fasta file for the mitochondria of all iberian females. 


```bash
# alignment with mafft
mafft --auto mitochondria.iberian_females.fa > mitochondria.iberian_females.mafft.fa

# triming alignments with trimal
trimal -in mitochondria.iberian_females.mafft.fa -out mitochondria.iberian_females.mafft.trimal.fa -gappyout
```

### Maximum likelihood tree construction with IQTREE-2

```bash
iqtree2 -s mitochondria.iberian_females.mafft.trimal.fa -m TEST -B 1000 -T AUTO --prefix mitochondrial_tree
```

## Neo-W 

The same pipeline was run with the neo-W VCF 

## Topology test

Topology testing with IQTREE-2 requires trimmed alignments, as well as the three whose topology needs to be tested (with trees the Newick format in columns).

```bash
iqtree2 -s mitochondria.iberian_females.mafft.trimal.fa -z topology_test.trees -zb 10000 -au -n 0 --prefix topology_test
```


