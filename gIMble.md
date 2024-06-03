# Demographic history inference

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

## MSMC

The demographic history of the neo-Z chromosome was inferred with MSMC2 in males 

```bash
#!/bin/bash

# msmc_pipeline.sh

ind=$1

bcftools view -s $ind -Oz -o vcfs/neosex.$ind.vcf.gz melanargia_ines.neosex.vcf.gz

# the multicallable file format is explained in the diversity_divergence notebook
grep $ind neosex.multicallable.bed > beds/$ind.bed

../msmc-tools/generate_multihetsep.py vcfs/neosex.$ind.vcf.gz --mask beds/$ind.bed > inputs/$ind.multi_het_sep.txt

# r of 3.4 was chosen based on the mutation and recombination rate of Heliconius melpomene
../msmc2/build/release/msmc2 -i 30 -o results/$ind.neosex -r 3.4 inputs/$ind.multi_het_sep.txt
```

```bash
while read ind; do bash msmc_pipeline.sh $ind; done < individuals.txt
```

The same commands were applied to autosomal data for all individuals in the dataset.
