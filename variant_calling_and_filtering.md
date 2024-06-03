# Variant callling and filtering

## Fastp

```bash
while read id; do fastp -p --detect_adapter_for_pe -i $id.1.fq.gz -I $id.2.fq.gz -o $id.trim.1.fastq.gz -O $id.trim.2.fastq.gz --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 --json $id.json --html $id.html --thread 10; done < individuals.txt
```

## Whole-genome alignment with BWA-MEM


```bash
while read ind; do bash bwa mem -t 16 -R '@RG\tID:$ind\tSM:$ind' 0_ref/melanargia_ines.PT_MI_8.v2_0.sequences.fasta 1_trimmed_reads/$ind.1.trim.fq.gz 1_trimmed_reads/$ind.2.trim.fq.gz | samtools view -q 1 -b - | samtools sort -@16 -m2G -T /scratch/tdecroly/$ind.tmp.bam - > 2_bwa_mem/$ind.vs.melanargia_ines.PT_MI_8.v2_0.sequences.bam && samtools index 2_bwa_mem/$ind.vs.melanargia_ines.PT_MI_8.v2_0.sequences.bam $ind; done < individuals.txt

# melanargia_ines.PT_MI_8.v2_0.sequences.fasta is the reference assembly
# individuals.txt is a text file containing individuals IDs in a column
```

## Marking PCR duplicates with Sambamba 

```bash
while read ind; do bash sambamba markdup -t 16 --tmpdir /scratch/tdecroly 2_bwa_mem/$ind.vs.melanargia_ines.PT_MI_8.v2_0.sequences.bam 2_bwa_mem/$ind.vs.melanargia_ines.PT_MI_8.v2_0.sequences.MD.bam $ind; done < individuals.txt
```

## Calling variants with freebayes 

```bash
freebayes-parallel <(fasta_generate_regions.py 0_ref/melanargia_ines.PT_MI_8.v2_0.sequences.fasta 100000000) 14 -f 0_ref/melanargia_ines.PT_MI_8.v2_0.sequences.fasta --limit-coverage 500 --use-best-n-alleles 8 --no-population-priors --use-mapping-quality --ploidy 2 --haplotype-length -1 --bam 2_bwa_mem/ES_MI_1647.vs.melanargia_ines.PT_MI_8.v2_0.sequences.MD.bam --bam 2_bwa_mem/PT_MI_7.vs.melanargia_ines.PT_MI_8.v2_0.sequences.MD.bam --bam 2_bwa_mem/PT_MI_8.vs.melanargia_ines.PT_MI_8.v2_0.sequences.MD.bam --bam 2_bwa_mem/PT_MI_61.vs.melanargia_ines.PT_MI_8.v2_0.sequences.MD.bam --bam 2_bwa_mem/MA_MI_1620.vs.melanargia_ines.PT_MI_8.v2_0.sequences.MD.bam --bam 2_bwa_mem/DZ_MI_1624.vs.melanargia_ines.PT_MI_8.v2_0.sequences.MD.bam --bam 2_bwa_mem/TN_MI_1619.vs.melanargia_ines.PT_MI_8.v2_0.sequences.MD.bam --bam 2_bwa_mem/PT_MI_86.vs.melanargia_ines.PT_MI_8.v2_0.sequences.MD.bam --bam 2_bwa_mem/ES_MI_1686.vs.melanargia_ines.PT_MI_8.v2_0.sequences.MD.bam --bam 2_bwa_mem/ES_MI_1684.vs.melanargia_ines.PT_MI_8.v2_0.sequences.MD.bam --bam 2_bwa_mem/ES_MI_1683.vs.melanargia_ines.PT_MI_8.v2_0.sequences.MD.bam --bam 2_bwa_mem/ES_MI_1682.vs.melanargia_ines.PT_MI_8.v2_0.sequences.MD.bam --bam 2_bwa_mem/ES_MI_1680.vs.melanargia_ines.PT_MI_8.v2_0.sequences.MD.bam --bam 2_bwa_mem/DZ_MI_1685.vs.melanargia_ines.PT_MI_8.v2_0.sequences.MD.bam --bam 2_bwa_mem/DZ_MI_1681.vs.melanargia_ines.PT_MI_8.v2_0.sequences.MD.bam > 3_freebayes/melanargia_ines.vs.melanargia_ines.PT_MI_8.variants_14.vcf.gz && bgzip 3_freebayes/melanargia_ines.vs.melanargia_ines.PT_MI_8.variants_14.vcf.gz
```

## Filtering variants with gIMble preprocess 

### Nuclear genome

```bash
gIMble preprocess -f 0_ref/melanargia_ines.PT_MI_8.v2_0.sequences.fasta -v 3_freebayes/melanargia_ines.vs.melanargia_ines.PT_MI_8.variants_14.vcf.gz -b /data/lohse/modelgenomland/analyses/melanargia_ines/2_bwa_mem -g 2 -q 10 -m 8 -M 3 -t 20 -o 4_preprocess/melanargia_ines.gimbleprepped -k
```

#### Indels-containing vcf

The indels-containing vcf the neo-sex chromosome was created with a custom commands largely inspired from gIMbleprep (cf [gIMble repo](https://github.com/DRL/gimbleprep)). The vcf was filter with:

```bash
bcftools norm neoZ.freebayes.vcf.gz -f ../fasta/melanargia_ines.PT_MI_8.v2_0.neoZ.fasta | vcfallelicprimitives --keep-info --keep-geno -t decomposed | bcftools filter -Ov -S . -e '(FMT/DP[0]<8 | FMT/DP[0]>=72) | (FMT/DP[1]<8 | FMT/DP[1]>=66) | (FMT/DP[2]<8 | FMT/DP[2]>=69) | (FMT/DP[3]<8 | FMT/DP[3]>=66) | (FMT/DP[4]<8 | FMT/DP[4]>=51) | (FMT/DP[5]<8 | FMT/DP[5]>=279) | (FMT/DP[6]<8 | FMT/DP[6]>=57) | (FMT/DP[7]<8 | FMT/DP[7]>=66) | (FMT/DP[8]<8 | FMT/DP[8]>=159) | (FMT/DP[9]<8 | FMT/DP[9]>=72) | (FMT/DP[10]<8 | FMT/DP[10]>=60) | (FMT/DP[11]<8 | FMT/DP[11]>=57) | (FMT/DP[12]<8 | FMT/DP[12]>=87) | (FMT/DP[13]<8 | FMT/DP[13]>=72) | (FMT/DP[14]<8 | FMT/DP[14]>=69)' | bcftools filter -Oz -i 'QUAL >= 10 & RPL>=1 & RPR>=1 & SAF>=1 & SAR>=1' | bcftools filter -Oz --SnpGap 2 -Oz -o melanargia_ines.neoZ.snps_structural.vcf.gz
```

Additionally, the multi-callable bed file was generated with the following commands: 

```bash
zgrep "melanargia_ines.PT_MI_8.chromosome_14" multi_callable.bed.gz | gzip > neoZ.multi_callable.bed.gz

bcftools norm neoZ.freebayes.vcf.gz -f ../fasta/melanargia_ines.PT_MI_8.v2_0.neoZ.fasta | vcfallelicprimitives --keep-info --keep-geno -t decomposed | bcftools filter -Oz -s Balance -m+ -i 'QUAL >= 10 & RPL>=1 & RPR>=1 & SAF>=1 & SAR>=1' | bcftools filter -Oz -m+ -s+ --SnpGap 2 -Oz -o tmp.neosex.structural_snp.filtering.vcf.gz

bcftools view -H -i "%FILTER!='PASS'" tmp.neosex.structural_snp.filtering.vcf.gz | perl -lane '$pad=0; print($F[0]."\t".($F[1]-1)."\t".(($F[1]-1)+length($F[3]))."\t".$F[6])' | bedtools sort | bedtools merge -i - | gzip > structural_snps.bed_fail.bed.gz

bedtools subtract -a ../bed/neoZ.multi_callable.bed.gz -b structural_snps.bed_fail.bed.gz | bedtools sort | gzip > neoZ.structural_snps.multi_callable.bed.gz
```

The `multi_callable.bed.gz` is found in the gIMble temporary directory. It contains region of the genome that falls within the coverage boundaries defined by the user. It was created by merging the `mosdepth --quantize` output of all the samples (that is, region with read depth between X and Y), using the `bedtools multiinter` command.  


## Read depth with mosdepth 

### Per base

```bash
while read ind; do mosdepth -t 4 -m 5_depth/$ind.vs.melanargia_ines.PT_MI_8.v2_0.sequences 2_bwa_mem/$ind.vs.melanargia_ines.PT_MI_8.v2_0.sequences.MD.bam ; done < individual.txt
```

### Per chromosome 

```bash
while read ind; do mosdepth -t 4 -m -b 99999999 -n 5_depth/$ind.vs.melanargia_ines.PT_MI_8.v2_0.sequences 2_bwa_mem/$ind.vs.melanargia_ines.PT_MI_8.v2_0.sequences.MD.bam; done < individual.txt
```

