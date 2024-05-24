I was generating genotype likelihoods earlier using angsd but apparently that sample list contained much more samples like *D. melanogaster*, *D. innubila* etc. and not just my *D. affinis* samples.

Here's the link for my [data ID file](https://github.com/anjaligupta1210/AG_Unckless_Open_Lab_Notebook/blob/master/datafiles/affinis_popgen_samples%20-%20Sheet1.pdf) and the [file with sex-ratio info](https://github.com/anjaligupta1210/AG_Unckless_Open_Lab_Notebook/blob/master/datafiles/affinis_pop_samples%20-SexRatio.csv).

I created a file called `bams_list.txt` containing my *D. affinis* pop gen files list -


```python
/work/unckless/a948g501/PopGen/Daff_ST.bam
/work/unckless/a948g501/PopGen/Daff_SR1.bam
/work/unckless/a948g501/PopGen/Daff_SR2.bam
/work/unckless/a948g501/PopGen/pool10.bam
/work/unckless/a948g501/PopGen/pool11.bam
/work/unckless/a948g501/PopGen/pool12.bam
/work/unckless/a948g501/PopGen/pool13.bam
/work/unckless/a948g501/PopGen/pool14.bam
/work/unckless/a948g501/PopGen/pool15.bam
/work/unckless/a948g501/PopGen/pool16.bam
/work/unckless/a948g501/PopGen/pool17.bam
/work/unckless/a948g501/PopGen/pool18.bam
/work/unckless/a948g501/PopGen/pool19.bam
/work/unckless/a948g501/PopGen/pool1.bam
/work/unckless/a948g501/PopGen/pool20.bam
/work/unckless/a948g501/PopGen/pool21.bam
/work/unckless/a948g501/PopGen/pool22.bam
/work/unckless/a948g501/PopGen/pool23.bam
/work/unckless/a948g501/PopGen/pool24.bam
/work/unckless/a948g501/PopGen/pool26.bam
/work/unckless/a948g501/PopGen/pool27.bam
/work/unckless/a948g501/PopGen/pool28.bam
/work/unckless/a948g501/PopGen/pool29.bam
/work/unckless/a948g501/PopGen/pool2.bam
/work/unckless/a948g501/PopGen/pool30.bam
/work/unckless/a948g501/PopGen/pool31.bam
/work/unckless/a948g501/PopGen/pool32.bam
/work/unckless/a948g501/PopGen/pool33.bam
/work/unckless/a948g501/PopGen/pool34.bam
/work/unckless/a948g501/PopGen/pool35.bam
/work/unckless/a948g501/PopGen/pool36.bam
/work/unckless/a948g501/PopGen/pool37.bam
/work/unckless/a948g501/PopGen/pool38.bam
/work/unckless/a948g501/PopGen/pool39.bam
/work/unckless/a948g501/PopGen/pool3.bam
/work/unckless/a948g501/PopGen/pool40.bam
/work/unckless/a948g501/PopGen/pool41.bam
/work/unckless/a948g501/PopGen/pool42.bam
/work/unckless/a948g501/PopGen/pool43.bam
/work/unckless/a948g501/PopGen/pool44.bam
/work/unckless/a948g501/PopGen/pool45.bam
/work/unckless/a948g501/PopGen/pool46.bam
/work/unckless/a948g501/PopGen/pool47.bam
/work/unckless/a948g501/PopGen/pool48.bam
/work/unckless/a948g501/PopGen/pool4.bam
/work/unckless/a948g501/PopGen/pool57.bam
/work/unckless/a948g501/PopGen/pool58.bam
/work/unckless/a948g501/PopGen/pool59.bam
/work/unckless/a948g501/PopGen/pool5.bam
/work/unckless/a948g501/PopGen/pool60.bam
/work/unckless/a948g501/PopGen/pool61.bam
/work/unckless/a948g501/PopGen/pool62.bam
/work/unckless/a948g501/PopGen/pool63.bam
/work/unckless/a948g501/PopGen/pool64.bam
/work/unckless/a948g501/PopGen/pool65.bam
/work/unckless/a948g501/PopGen/pool66.bam
/work/unckless/a948g501/PopGen/pool67.bam
/work/unckless/a948g501/PopGen/pool68.bam
/work/unckless/a948g501/PopGen/pool69.bam
/work/unckless/a948g501/PopGen/pool6.bam
/work/unckless/a948g501/PopGen/pool70.bam
/work/unckless/a948g501/PopGen/pool71.bam
/work/unckless/a948g501/PopGen/pool72.bam
/work/unckless/a948g501/PopGen/pool74.bam
/work/unckless/a948g501/PopGen/pool75.bam
/work/unckless/a948g501/PopGen/pool76.bam
/work/unckless/a948g501/PopGen/pool77.bam
/work/unckless/a948g501/PopGen/pool78.bam
/work/unckless/a948g501/PopGen/pool79.bam
/work/unckless/a948g501/PopGen/pool7.bam
/work/unckless/a948g501/PopGen/pool80.bam
/work/unckless/a948g501/PopGen/pool87.bam
/work/unckless/a948g501/PopGen/pool88.bam
/work/unckless/a948g501/PopGen/pool89.bam
/work/unckless/a948g501/PopGen/pool8.bam
/work/unckless/a948g501/PopGen/pool90.bam
/work/unckless/a948g501/PopGen/pool91.bam
/work/unckless/a948g501/PopGen/pool9.bam

```

### Generate genotype likelihoods and pcangsd on ChrX

Next I generated genotype likelihoods on my *D. affinis* samples for ChrX using angsd


```python
#!/bin/bash
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=100:00:00
#SBATCH --mem=100GB
#SBATCH --partition=unckless
#SBATCH --chdir=/work/unckless/a948g501/PopGen/    
#SBATCH -o sh.%j.out
#SBATCH -e sh.%j.err

module load angsd


# Generate genotype likelihoods
angsd -bam bam_list.txt -r ChrX_MullerAD -ref Daffinis_STfemale_v5.1.masked.fasta -GL 2 -doGlf 2 -doMajorMinor 1 -doCounts 1 -doHaploCall 2 -doMaf 1 -minMaf 0.05 -out output_prefix_X
```

PCangsd on genotype likelihoods for ChrX -


```python
#!/bin/bash
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=06:00:00
#SBATCH --mem=100GB
#SBATCH --partition=sixhour
#SBATCH -o sh.%j.out
#SBATCH -e sh.%j.err

module load python
pcangsd -b output_prefix_X.beagle.gz -e 2 -t 64 -o pcangsd_X --selection
```

### Generate genotype likelihoods and pcangsd on Autosomes

Next I generated genotype likelihoods on my *D. affinis* samples for autosomes using angsd


```python
#!/bin/bash
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=100:00:00
#SBATCH --mem=100GB
#SBATCH --partition=unckless
#SBATCH --chdir=/work/unckless/a948g501/PopGen/    
#SBATCH -o sh.%j.out
#SBATCH -e sh.%j.err


# Specify the bed file containing regions to include
bed_file="Autosomes.bed"

module load angsd

angsd -bam bam_list.txt -rf $bed_file -ref Daffinis_STfemale_v5.1.masked.fasta -GL 2 -doGlf 2 -doMajorMinor 1 -doCounts 1 -doMaf 1 -minMaf 0.05 -out output_prefix_Autosomes
```

PCangsd on genotype likelihoods for autosomes -


```python
#!/bin/bash
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=06:00:00
#SBATCH --mem=500GB
#SBATCH --partition=sixhour
#SBATCH -o sh.%j.out
#SBATCH -e sh.%j.err

module load python
pcangsd -b output_prefix_Autosomes.beagle.gz -e 2 -t 64 -o pcangsd_Autosome --selection
```

### Generate genotype likelihoods and pcangsd separately for autosomes - Chr2, 3, 4

**Chr 4 -**


```python
#!/bin/bash
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=100:00:00
#SBATCH --mem=100GB
#SBATCH --partition=unckless
#SBATCH --chdir=/work/unckless/a948g501/PopGen/    
#SBATCH -o sh.%j.out
#SBATCH -e sh.%j.err




module load angsd

angsd -bam bam_list.txt -r Chr4_MullerB -ref Daffinis_STfemale_v5.1.masked.fasta -GL 2 -doGlf 2 -doMajorMinor 1 -doCounts 1 -doMaf 1 -minMaf 0.05 -out output_prefix_Chr4_MullerB

module load python

pcangsd -b output_prefix_Chr4_MullerB.beagle.gz -e 2 -t 64 -o pcangsd_Chr4_MullerB --selection
```

**Chr 3 -**


```python
#!/bin/bash
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=100:00:00
#SBATCH --mem=100GB
#SBATCH --partition=unckless
#SBATCH --chdir=/work/unckless/a948g501/PopGen/    
#SBATCH -o sh.%j.out
#SBATCH -e sh.%j.err




module load angsd

angsd -bam bam_list.txt -r Chr3_MullerC -ref Daffinis_STfemale_v5.1.masked.fasta -GL 2 -doGlf 2 -doMajorMinor 1 -doCounts 1 -doMaf 1 -minMaf 0.05 -out output_prefix_Chr3_MullerC

module load python

pcangsd -b output_prefix_Chr3_MullerC.beagle.gz -e 2 -t 64 -o pcangsd_Chr3_MullerC --selection
```

**Chr 2-**

Chr2.bed file -


```python
Chr2_MullerE    1    44456471
Chr2.group4_MullerE    1    50000
```


```python
#!/bin/bash
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=100:00:00
#SBATCH --mem=100GB
#SBATCH --partition=unckless
#SBATCH --chdir=/work/unckless/a948g501/PopGen/    
#SBATCH -o sh.%j.out
#SBATCH -e sh.%j.err


# Specify the bed file containing regions to include
bed_file="Chr2.bed"

module load angsd

angsd -bam bam_list.txt -rf $bed_file -ref Daffinis_STfemale_v5.1.masked.fasta -GL 2 -doGlf 2 -doMajorMinor 1 -doCounts 1 -doMaf 1 -minMaf 0.05 -out output_prefix_Chr2

module load python

pcangsd -b output_prefix_Chr2.beagle.gz -e 2 -t 64 -o pcangsd_Chr2 --selection
```
