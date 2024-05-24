## **Calling SNPs using GATK**

To call SNPs on my *D. affinis* ST, SR1, SR2 and other related species with only one sample per species, so far I used -

1. bcftools mpileup - Not great because the pipeline is not sensitive enough to capture all variation, so we didn't know if our vcf is missing important variation across species.
2. Freebayes - Pipeline is very sensitive and captures all variation but the vcf files produced by freebayes are weird because when I try to convert them to phylip alignments, it somehow only shows up to be a single sample only no matter how many samples I have in my vcf's

So now I'm moving to GATK to call SNPs, Fingers crossed this works!

Using GATK -

First, I need to index my reference fasta file with a .dict extension using gatk CreateSequenceDictionary


```python
gatk CreateSequenceDictionary -R Daffinis_STfemale_v5.1.masked.fasta
```

Next, for each bam alignment file, we need to use Picard to AddOrReplaceReadGroups before we can use it to call variants

Here's my script to do so, it is called PicardAddOrReplaceReadGroups.sh


```python
#!/bin/bash
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=06:00:00
#SBATCH --mem=10GB
#SBATCH --partition=sixhour
#SBATCH -o sh.%j.out
#SBATCH -e sh.%j.err

i=$1
module load picard
picard AddOrReplaceReadGroups I=${i}.bam O=${i}.picard.bam RGLB=${i} RGID=${i} RGPL=illumina RGPU=unit1 RGSM=${i}

```

I'm just going to run these in a loop so that it is faster and each bam alignment gets submitted as a separate job


```python
for i in "Daff_SR1" "Daff_SR2" "Daff_ST" "Dalg" "Datha_ea" "Datha_eb" "Dazt" "Dpse"
do
sbatch PicardAddOrReplaceReadGroups.sh $i
done
```

Now, we are going to index our post-PicardAddOrReplaceReadGroups bam alignment files using samtools


```python
module load samtools

for i in "Daff_SR1" "Daff_SR2" "Daff_ST" "Dalg" "Datha_ea" "Datha_eb" "Dazt" "Dpse"
do
nohup samtools index ${i}.picard.bam &
done
```

Now, I'm going to use gatk HaplotypeCaller to call SNPs separately for ChrX (haploid) and autosomes (diploid)

Here are the two scripts to call SNPs using gatk HaplotypeCaller


quick things: 

- -L is to specify only region to include
- -XL is to specify only region to exclude


GATK is slow and will take a long time to run 


```python
#!/bin/bash
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=60-00:00:00
#SBATCH --mem=100GB
#SBATCH --partition=unckless
#SBATCH -o sh.%j.out
#SBATCH -e sh.%j.err


module load gatk
gatk HaplotypeCaller \
     -I Daff_SR1.picard.bam \
     -I Daff_SR2.picard.bam \
     -I Daff_ST.picard.bam \
     -I Dalg.picard.bam \
     -I Datha_ea.picard.bam \
     -I Datha_eb.picard.bam \
     -I Dazt.picard.bam \
     -I Dpse.picard.bam \
     -O GATK_ChrX_haploid.vcf \
     -R Daffinis_STfemale_v5.1.masked.fasta \
     --add-output-vcf-command-line true \
     -ploidy 1 \
     -L ChrX_MullerAD

```


```python
#!/bin/bash
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=60-00:00:00
#SBATCH --mem=100GB
#SBATCH --partition=unckless
#SBATCH -o sh.%j.out
#SBATCH -e sh.%j.err



module load gatk
gatk HaplotypeCaller \
     -I Daff_SR1.picard.bam \
     -I Daff_SR2.picard.bam \
     -I Daff_ST.picard.bam \
     -I Dalg.picard.bam \
     -I Datha_ea.picard.bam \
     -I Datha_eb.picard.bam \
     -I Dazt.picard.bam \
     -I Dpse.picard.bam \
     -O GATK_Autosomes_diploid.vcf \
     -R Daffinis_STfemale_v5.1.masked.fasta \
     --add-output-vcf-command-line true \
     -ploidy 2 \
     -XL ChrX_MullerAD
```

Next I wanted to do some QC on my vcf's from gatk to filter only high quality varriants/ SNPs, so I used bcftools filter -


```python
module load bcftools
bcftools filter -i 'QUAL > 40' GATK_ChrX_haploid.vcf -o GATK_filtered_ChrX_haploid.vcf

bcftools filter -i 'QUAL > 40 && TYPE="snp"' GATK_ChrX_haploid.vcf -o GATK_filtered_ChrX_snp_haploid.vcf
```


```python
module load bcftools
bcftools filter -i 'QUAL > 40 && TYPE="snp"' GATK_Autosomes_diploid.vcf -o GATK_filtered_Autosome_snp_diploid.vcf
```

## **Making Trees for Chromosome X**

First I converted my all variants filtered ChrX vcf to a phylip alignment -


```python
module load python 
python vcf2phylip.py -i GATK_filtered_ChrX_haploid.vcf
```

### Making trees across ChrX in sliding windows of 1 Kb

GATK_filtered_ChrX_haploid.min4.phy

Length of the phylip alignment is 1,935,737

So I will generate 1935 1Kb alignments and ignore the last 737 bases

To do so, I will define two variables in a for loop.

One (i) will track 1-1935, over the 1935 iterations one (j) will track 1-1,934,001 in increments of 1Kb, over the 1935 iterations

The code ‘$((j+999))’ means add 999 to j each iteration, which will track 1,000-1,935,000 in increments of 1Kb, over the 1935 iterations

Run this code to generate the necessary partition file


```python
j=1
for i in {1..1935}
do
echo "DNA, p${i}=${j}-$((j+999))" >> partitions_X_GATK.txt
j=$((j+1000))
done
```


```python
mkdir partitioned_X_GATK.phylips
module load python
python unconcatenate_phylip.py GATK_filtered_ChrX_haploid.min4.phy partitions_X_GATK.txt --prefix=partitioned_X_GATK.phylips/
```


```python
#!/bin/sh
#
#SBATCH --job-name=auto.gene.trees.X               
#SBATCH --nodes=1             
#SBATCH --ntasks-per-node=15               
#SBATCH --partition=sixhour            
#SBATCH --chdir=/work/unckless/a948g501/SlidingTrees/partitioned_X_GATK.phylips    
#SBATCH --mem-per-cpu=2gb            
#SBATCH --array=1-1935
#SBATCH --time=360


#To be able to later use bootstrapping in ASTRAL, we will also generate sets of trees from bootstrapped gene alignments in the same IQ-TREE analysis, by using the option -B 
#we do not specify a substition model with the -m option and therefore allow IQ-TREE to automatically select the best-fitting model.
#we will now ensure that the bootstrap trees are actually written to a file, and not just summarized in the output, by specifying the option --wbt.
#-nt AUTO ensures that all available CPU's are used

##load module iqtree
module load iqtree
iqtree2 -s /work/unckless/a948g501/SlidingTrees/partitioned_X_GATK.phylips/DNA_p$SLURM_ARRAY_TASK_ID.phylip -bb 1000 -wbt -nt AUTO
```


```python
cat *.treefile > outputTrees_GATK_X.trees
```

### Making Trees across ChrX in sliding windows of 10kb

GATK_filtered_ChrX_haploid.min4.phy

Length of the phylip alignment is 1,935,737

So I will generate 193 10Kb alignments and ignore the last 5,737 bases

To do so, I will define two variables in a for loop.

One (i) will track 1-193, over the 193 iterations one (j) will track 1-1,920,001 in increments of 10Kb, over the 193 iterations

The code ‘$((j+9999))’ means add 9999 to j each iteration, which will track 10,000-1,930,000 in increments of 10Kb, over the 193 iterations

Run this code to generate the necessary partition file


```python
j=1
for i in {1..193}
do
echo "DNA, p${i}=${j}-$((j+9999))" >> partitions_X_GATK_10kb.txt
j=$((j+10000))
done
```


```python
mkdir partitioned_X_GATK_10kb.phylips
module load python
python unconcatenate_phylip.py GATK_filtered_ChrX_haploid.min4.phy partitions_X_GATK_10kb.txt --prefix=partitioned_X_GATK_10kb.phylips/
```


```python
#!/bin/sh
#
#SBATCH --job-name=auto.gene.trees.X               
#SBATCH --nodes=1             
#SBATCH --ntasks-per-node=15               
#SBATCH --partition=sixhour            
#SBATCH --chdir=/work/unckless/a948g501/SlidingTrees/partitioned_X_GATK_10kb.phylips    
#SBATCH --mem-per-cpu=2gb            
#SBATCH --array=1-193
#SBATCH --time=360


#To be able to later use bootstrapping in ASTRAL, we will also generate sets of trees from bootstrapped gene alignments in the same IQ-TREE analysis, by using the option -B 
#we do not specify a substition model with the -m option and therefore allow IQ-TREE to automatically select the best-fitting model.
#we will now ensure that the bootstrap trees are actually written to a file, and not just summarized in the output, by specifying the option --wbt.
#-nt AUTO ensures that all available CPU's are used

##load module iqtree
module load iqtree
iqtree2 -s /work/unckless/a948g501/SlidingTrees/partitioned_X_GATK_10kb.phylips/DNA_p$SLURM_ARRAY_TASK_ID.phylip -bb 1000 -wbt -nt AUTO
```


```python
cat *.treefile > outputTrees_GATK_X_10kb.trees
```

### Making trees across ChrX in sliding windows of 1 kb - but using only SNPs vcf instead of "all-variants" vcf


```python
module load python 
python vcf2phylip.py -i GATK_filtered_ChrX_snp_haploid.vcf
```

GATK_filtered_ChrX_haploid.min4.phy

Length of the phylip alignment is 1,935,737

So I will generate 1935 1Kb alignments and ignore the last 737 bases

To do so, I will define two variables in a for loop.

One (i) will track 1-1935, over the 1935 iterations one (j) will track 1-1,934,001 in increments of 1Kb, over the 1935 iterations

The code ‘$((j+999))’ means add 999 to j each iteration, which will track 1,000-1,935,000 in increments of 1Kb, over the 1935 iterations

Run this code to generate the necessary partition file


```python
j=1
for i in {1..1935}
do
echo "DNA, p${i}=${j}-$((j+999))" >> partitions_X_GATK.txt
j=$((j+1000))
done
```


```python
mkdir partitioned_X_GATK_snp.phylips
module load python
python unconcatenate_phylip.py GATK_filtered_ChrX_snp_haploid.min4.phy partitions_X_GATK.txt --prefix=partitioned_X_GATK_snp.phylips/
```


```python
#!/bin/sh
#
#SBATCH --job-name=auto.gene.trees.X               
#SBATCH --nodes=1             
#SBATCH --ntasks-per-node=15               
#SBATCH --partition=sixhour            
#SBATCH --chdir=/work/unckless/a948g501/SlidingTrees/partitioned_X_GATK_snp.phylips    
#SBATCH --mem-per-cpu=2gb            
#SBATCH --array=1-1935
#SBATCH --time=360


#To be able to later use bootstrapping in ASTRAL, we will also generate sets of trees from bootstrapped gene alignments in the same IQ-TREE analysis, by using the option -B 
#we do not specify a substition model with the -m option and therefore allow IQ-TREE to automatically select the best-fitting model.
#we will now ensure that the bootstrap trees are actually written to a file, and not just summarized in the output, by specifying the option --wbt.
#-nt AUTO ensures that all available CPU's are used

##load module iqtree
module load iqtree
iqtree2 -s /work/unckless/a948g501/SlidingTrees/partitioned_X_GATK_snp.phylips/DNA_p$SLURM_ARRAY_TASK_ID.phylip -bb 1000 -wbt -nt AUTO
```


```python
cat *.treefile > outputTrees_GATK_X_snp.trees
```

### Making trees across ChrX in sliding windows of 10kb - but using only SNPs vcf instead of "all-variants" vcf

GATK_filtered_ChrX_haploid.min4.phy

Length of the phylip alignment is 1,935,737

So I will generate 19 100Kb alignments and ignore the last 35,737 bases

To do so, I will define two variables in a for loop.

One (i) will track 1-19, over the 19 iterations one (j) will track 1-1,800,001 in increments of 100Kb, over the 19 iterations

The code ‘$((j+99999))’ means add 99,999 to j each iteration, which will track 1,00,000-1,900,000 in increments of 100Kb, over the 19 iterations

Run this code to generate the necessary partition file


```python
j=1
for i in {1..19}
do
echo "DNA, p${i}=${j}-$((j+99999))" >> partitions_X_GATK_snp_100kb.txt
j=$((j+100000))
done
```


```python
mkdir partitioned_X_GATK_snp_100kb.phylips
module load python
python unconcatenate_phylip.py GATK_filtered_ChrX_snp_haploid.min4.phy partitions_X_GATK_snp_100kb.txt --prefix=partitioned_X_GATK_snp_100kb.phylips/
```


```python
#!/bin/sh
#
#SBATCH --job-name=auto.gene.trees.X               
#SBATCH --nodes=1             
#SBATCH --ntasks-per-node=15               
#SBATCH --partition=sixhour            
#SBATCH --chdir=/work/unckless/a948g501/SlidingTrees/partitioned_X_GATK_snp_100kb.phylips    
#SBATCH --mem-per-cpu=2gb            
#SBATCH --array=1-19
#SBATCH --time=360


#To be able to later use bootstrapping in ASTRAL, we will also generate sets of trees from bootstrapped gene alignments in the same IQ-TREE analysis, by using the option -B 
#we do not specify a substition model with the -m option and therefore allow IQ-TREE to automatically select the best-fitting model.
#we will now ensure that the bootstrap trees are actually written to a file, and not just summarized in the output, by specifying the option --wbt.
#-nt AUTO ensures that all available CPU's are used

##load module iqtree
module load iqtree
iqtree2 -s /work/unckless/a948g501/SlidingTrees/partitioned_X_GATK_snp_100kb.phylips/DNA_p$SLURM_ARRAY_TASK_ID.phylip -bb 1000 -wbt -nt AUTO
```


```python
cat *.treefile > outputTrees_GATK_X_snp_100kb.trees
```

### Tree for the entire X


```python
iqtree2 -s GATK_filtered_ChrX_snp_haploid.min4.phy -bb 1000 -wbt -nt AUTO
```
