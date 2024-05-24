I made trees in sliding windows after calling SNPs using bcftools mpileup, but apparently bcftools mpileup is not the best tool to call SNPs when I only have one sample per species.
So I'm going to call SNPs on my bam files using freebayes and make trees in sliding windows again.

Here is my script to call SNPs on my alignments using freebayes - I'm calling it CallVariantsFreebayes.sh


```python
#!/bin/sh
#
#SBATCH --job-name=auto.freebayes              
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=60-00:00:00
#SBATCH --mem=100GB
#SBATCH --partition=unckless
#SBATCH -o sh.%j.out
#SBATCH -e sh.%j.err          
#SBATCH --chdir=/work/unckless/a948g501/SlidingTrees/  


# Load module
module load freebayes
module load bcftools

# Haploid calling for ChrX_MullerAD
freebayes -f Daffinis_STfemale_v5.1.masked.fasta --ploidy 1 --bam Daff_SR1.bam --bam Daff_SR2.bam --bam Daff_ST.bam --bam Dalg.bam --bam Datha_ea.bam --bam Datha_eb.bam --bam Dazt.bam --bam Dpse.bam -r ChrX_MullerAD > freebayes_ChrX_haploid.vcf

bcftools filter -i 'QUAL > 30' freebayes_ChrX_haploid.vcf -o filtered_freebayes_ChrX_haploid.vcf


# Diploid calling for autosomes
freebayes -f Daffinis_STfemale_v5.1.masked.fasta --ploidy 2 --bam Daff_SR1.bam --bam Daff_SR2.bam --bam Daff_ST.bam --bam Dalg.bam --bam Datha_ea.bam --bam Datha_eb.bam --bam Dazt.bam --bam Dpse.bam -r Chr4_MullerB -r Chr2_MullerE -r Chr2.group4_MullerE -r Chr3_MullerC -r Chr5_MullerF -r Unknown_69 -r mtDNA > freebayes_Autosomes_diploid.vcf

bcftools filter -i 'QUAL > 30' freebayes_Autosomes_diploid.vcf -o filtered_freebayes_Autosomes_diploid.vcf

```

I'll make my script executable and then run it using sbatch


```python
chmod +x CallVariantsFreebayes.sh
```

Now, I got the all-variant vcf files from Freebayes and I'll convert these into phylip alignments to make trees


```python
module load python 
python vcf2phylip.py -i filtered_freebayes_ChrX_haploid.vcf
python vcf2phylip.py -i filtered_freebayes_Autosomes_diploid.vcf
```

So when I was converting these vcf's from Freebayes into phylips, for some weird reason, it says number of samples in my vcf = 1 instead of 8 and I ran this thing twice and the vcf looks like it has all 8 samples but something is going wrong so I'm not going to use Freebayes to call SNPs.
