I want to call SNPs from my alignments following mapping the reads from 90 male samples to the D. affinis female genome.

I am using angsd to call SNPs from these alignments.

I'm doing this on the cluster in /work/unckless/a948g501/PopGen/

I will first generate genotype likelihoods, and then call and filter SNPs using angsd. At the end, I will convert it into a vcf format using angsd.

I will do this separately for the X chromosome and then for the autosomes. I will also use this pipeline later for the Y chromosome.

First, I generated a text file of my bam files list.

My pipeline to generate genotype likelihoods is called angsdGenerateGenotypeLikelihoods.sh and first I'm using this only for the X chromosome. Here is the script for angsdGenerateGenotypeLikelihoods.sh


```python
#!/bin/bash

module load angsd

ls /work/unckless/a948g501/PopGen/*.bam > bam_list.txt

# Generate genotype likelihoods
angsd -bam bam_list.txt -r ChrX_MullerAD -ref Daffinis_STfemale_v5.1.fasta -GL 2 -doGlf 2 -doMajorMinor 1 -doCounts 1 -doHaploCall 2 -out output_prefix

```

Next, I made my script executable using the following


```python
chmod +x angsdGenerateGenotypeLikelihoods.sh
```

Then, I ran the script on the cluster using the following Script.sh


```python
#!/bin/bash
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=100:00:00
#SBATCH --mem=100GB
#SBATCH --partition=unckless
#SBATCH -o sh.%j.out
#SBATCH -e sh.%j.err


sh angsdGenerateGenotypeLikelihoods.sh
```
