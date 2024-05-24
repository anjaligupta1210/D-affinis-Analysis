Rob suggested I do something different instead of doing regular Sliding window trees

So, what I'm doing now is 
1. tree for the entire X 
2. left arm vs right arm 
3. left and right arm partitioned using inversion breakpoints
4. Sliding window trees but instead of partititioning the phylip - which happens randomly - I'm partitioning the vcf file and converting separate vcfs into phylips into trees

We will use only SNP filtered files - GATK_filtered_ChrX_snp_haploid.vcf and GATK_filtered_Autosome_snp_diploid.vcf

Our azteca sample from NCBI is grouping weirdly for the entire autosome or entire X tree so we will exclude azteca from our vcf's for now!

### Getting rid of *D. azteca* from our vcf

We will get rid of azteca and blanks from our vcf file 

We will use bcftools view to remove azteca and SNP's from our vcf where something is missing in any one sample


```python
module load bcftools

bcftools view -s ^Dazt -e 'FMT/DP="." | FMT/AD="." | FMT/GQ="."' -o GATKfilteredSNPchrX_NoAztNoBlanks.vcf GATK_filtered_ChrX_snp_haploid.vcf
```


```python
module load bcftools

bcftools view -s ^Dazt -e 'FMT/DP="." | FMT/AD="." | FMT/GQ="."' -o GATKfilteredSNPAutosome_NoAztNoBlanks.vcf GATK_filtered_Autosome_snp_diploid.vcf
```

### Tree for all autosomes

Tree for all Autosomes - first convert vcf to phylip using [vcf2phylip.py](https://raw.githubusercontent.com/edgardomortiz/vcf2phylip/master/vcf2phylip.py) and then I'll make the tree using IQ-TREE2


```python
module load python
python vcf2phylip.py -i GATKfilteredSNPAutosome_NoAztNoBlanks.vcf 

module load iqtree
iqtree2 -s GATKfilteredSNPAutosome_NoAztNoBlanks.min4.phy -bb 1000 -wbt -nt AUTO
```

### Tree for the entire X

Tree for entire X - first convert vcf to phylip using [vcf2phylip.py](https://raw.githubusercontent.com/edgardomortiz/vcf2phylip/master/vcf2phylip.py) and then I'll make the tree using IQ-TREE2


```python
module load python
python vcf2phylip.py -i GATKfilteredSNPchrX_NoAztNoBlanks.vcf

module load iqtree
iqtree2 -s GATKfilteredSNPchrX_NoAztNoBlanks.min4.phy -bb 1000 -wbt -nt AUTO
```

### Tree for left arm and right arm of ChrX

Partition vcf for ChrX left and right arm - 

First I need to zip and index my vcf using samtools bgzip and tabix


```python
module load samtools
bgzip GATKfilteredSNPchrX_NoAztNoBlanks.vcf
tabix GATKfilteredSNPchrX_NoAztNoBlanks.vcf.gz
```

Now I will partition my vcf file into left and right arm sections assuming the centromere is at 30 Mb


```python
module load bcftools
bcftools view -r ChrX_MullerAD:1-30000000 GATKfilteredSNPchrX_NoAztNoBlanks.vcf.gz -o chrX_XL.vcf
bcftools view -r ChrX_MullerAD:30000000- GATKfilteredSNPchrX_NoAztNoBlanks.vcf.gz -o chrX_XR.vcf
```

Convert vcf to phylip using [vcf2phylip.py](https://raw.githubusercontent.com/edgardomortiz/vcf2phylip/master/vcf2phylip.py) and then I'll make the tree using IQ-TREE2

module load python
python vcf2phylip.py -i chrX_XL.vcf
python vcf2phylip.py -i chrX_XR.vcf

module load iqtree
iqtree2 -s chrX_XL.min4.phy -bb 1000 -wbt -nt AUTO
iqtree2 -s chrX_XR.min4.phy -bb 1000 -wbt -nt AUTO


### Trees for left and right arm based on inversion breakpoint positions

Partitioning vcf for the two single inversions using inversion breakpoint positions that Rob gave me


```python
module load bcftools
bcftools view -r ChrX_MullerAD:1-5490855 GATKfilteredSNPchrX_NoAztNoBlanks.vcf.gz -o chrX_XL_PreInv.vcf
bcftools view -r ChrX_MullerAD:5490855-10233152 GATKfilteredSNPchrX_NoAztNoBlanks.vcf.gz -o chrX_XL_Inv.vcf
bcftools view -r ChrX_MullerAD:10233152-30000000 GATKfilteredSNPchrX_NoAztNoBlanks.vcf.gz -o chrX_XL_PostInv.vcf

bcftools view -r ChrX_MullerAD:30000000-44023827 GATKfilteredSNPchrX_NoAztNoBlanks.vcf.gz -o chrX_XR_PreInv.vcf
bcftools view -r ChrX_MullerAD:44023827-62734159 GATKfilteredSNPchrX_NoAztNoBlanks.vcf.gz -o chrX_XR_Inv.vcf
bcftools view -r ChrX_MullerAD:62734159- GATKfilteredSNPchrX_NoAztNoBlanks.vcf.gz -o chrX_XR_PostInv.vcf
```

Convert vcf to phylip using [vcf2phylip.py](https://raw.githubusercontent.com/edgardomortiz/vcf2phylip/master/vcf2phylip.py) and then I'll make the tree using IQ-TREE2


```python
module load python
python vcf2phylip.py -i chrX_XL_PreInv.vcf
python vcf2phylip.py -i chrX_XL_Inv.vcf
python vcf2phylip.py -i chrX_XL_PostInv.vcf

python vcf2phylip.py -i chrX_XR_PreInv.vcf
python vcf2phylip.py -i chrX_XR_Inv.vcf
python vcf2phylip.py -i chrX_XR_PostInv.vcf

module load iqtree
iqtree2 -s chrX_XL_PreInv.min4.phy -bb 1000 -wbt -nt AUTO
iqtree2 -s chrX_XL_Inv.min4.phy -bb 1000 -wbt -nt AUTO
iqtree2 -s chrX_XL_PostInv.min4.phy -bb 1000 -wbt -nt AUTO

iqtree2 -s chrX_XR_PreInv.min4.phy -bb 1000 -wbt -nt AUTO
iqtree2 -s chrX_XR_Inv.min4.phy -bb 1000 -wbt -nt AUTO
iqtree2 -s chrX_XR_PostInv.min4.phy -bb 1000 -wbt -nt AUTO


```

### Trees in sliding windows of 100kb across the X (partitioning vcf's instead of the phylip)

Sliding trees across X

Partitioning my vcf files into multiple files 

First I'm going to make a bed file with information about each partition 


```python
j=1
for i in {1..680}
do
    echo -e "ChrX_MullerAD\t${j}\t$((j+99999))" >> partitions_X.bed
    j=$((j+100000))
done
```

Next I'll run bcftools view looping over each line in my partitions file


```python
#!/bin/bash

# Define the input BED file
bed_file="partitions_X.bed"

# Define the input VCF file
vcf_file="GATKfilteredSNPchrX_NoAztNoBlanks.vcf.gz"

module load bcftools

# Loop through each line of the BED file
while IFS=$'\t' read -r chromosome start end; do
    # Define the output file prefix
    output_prefix="${chromosome}_${start}_${end}"

    # Define the output VCF file name
    output_vcf="${output_prefix}.vcf"

    # Subset the VCF file using bcftools view
    bcftools view -r "${chromosome}:${start}-${end}" "$vcf_file" -o "$output_vcf"

    echo "Subset VCF file generated: $output_vcf"
done < "$bed_file"

```

Convert all vcfs to phylips


```python
#!/bin/bash

module load python

# Loop through each VCF file matching the pattern
for vcf_file in ChrX_MullerAD_*.vcf; do
    python vcf2phylip.py -i "$vcf_file"
done
```

Make trees for each alignment


```python
#!/bin/bash

module load iqtree

# Loop through each VCF file matching the pattern
for phy_file in ChrX_MullerAD_*.min4.phy; do
    iqtree2 -s "$phy_file" -bb 1000 -wbt -nt AUTO 
done
```
