---
layout: post
title: Trees in sliding windows of 1-4MB
date: 13 June 2024  
category: [ Computational Pipelines ]
tags: [ Sliding window trees, Comparative genomics ]
---

# Trees in sliding windows of 1-4MB across the X of *D. affinis* **ST,SR1,SR2** , *D. algonquin, D. athabasca, D. pseudoobscura*

## Trees in sliding windows of 1MB across the X (partitioning vcf’s instead of the phylip)

Sliding window trees across X

Partitioning my vcf files into multiple files

First I’m going to make a bed file with information about each partition


```python
j=1
for i in {1..68}
do
    echo -e "ChrX_MullerAD\t${j}\t$((j+999999))" >> partitions_X_1mb.bed
    j=$((j+1000000))
done
```

Next I’ll run bcftools view looping over each line in my partitions file


```python
#!/bin/bash

# Define the input BED file
bed_file="partitions_X_1mb.bed"

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

bed_file="partitions_X_1mb.bed"

module load python

# Loop through each VCF file matching the pattern

while IFS=$'\t' read -r chromosome start end; do
    # Define the output file prefix
    output_prefix="${chromosome}_${start}_${end}"

    for vcf_file in "${output_prefix}.vcf"; do
    python vcf2phylip.py -i "$vcf_file"
    done
done < "$bed_file"
```

Make trees for each alignment


```python
#!/bin/bash

bed_file="partitions_X_1mb.bed"

mkdir 1mbtrees_X

module load iqtree

# Loop through each VCF file matching the pattern

while IFS=$'\t' read -r chromosome start end; do
    # Define the output file prefix
    output_prefix="${chromosome}_${start}_${end}"

    # Loop through each VCF file matching the pattern
    for phy_file in "${output_prefix}.min4.phy"; do
    iqtree2 -s "$phy_file" -bb 1000 -wbt -nt AUTO 
    done

mv "$phy_file"* 1mbtrees_X

done < "$bed_file"
```

## Trees in sliding windows of 2MB across the X (partitioning vcf’s instead of the phylip)

Sliding window trees across X

Partitioning my vcf files into multiple files

First I’m going to make a bed file with information about each partition


```python
j=1
for i in {1..34}
do
    echo -e "ChrX_MullerAD\t${j}\t$((j+1999999))" >> partitions_X_2mb.bed
    j=$((j+2000000))
done
```

Next I’ll run bcftools view looping over each line in my partitions file


```python
#!/bin/bash

# Define the input BED file
bed_file="partitions_X_2mb.bed"

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

bed_file="partitions_X_2mb.bed"

module load python

# Loop through each VCF file matching the pattern

while IFS=$'\t' read -r chromosome start end; do
    # Define the output file prefix
    output_prefix="${chromosome}_${start}_${end}"

    for vcf_file in "${output_prefix}.vcf"; do
    python vcf2phylip.py -i "$vcf_file"
    done
done < "$bed_file"
```

Make trees for each alignment


```python
#!/bin/bash

bed_file="partitions_X_2mb.bed"

mkdir 2mbtrees_X

module load iqtree

# Loop through each VCF file matching the pattern

while IFS=$'\t' read -r chromosome start end; do
    # Define the output file prefix
    output_prefix="${chromosome}_${start}_${end}"

    # Loop through each VCF file matching the pattern
    for phy_file in "${output_prefix}.min4.phy"; do
    iqtree2 -s "$phy_file" -bb 1000 -wbt -nt AUTO 
    done

mv "$phy_file"* 2mbtrees_X

done < "$bed_file"
```

## Trees in sliding windows of 3MB across the X (partitioning vcf’s instead of the phylip)

Sliding window trees across X

Partitioning my vcf files into multiple files

First I’m going to make a bed file with information about each partition


```python
j=1
for i in {1..23}
do
    echo -e "ChrX_MullerAD\t${j}\t$((j+2999999))" >> partitions_X_3mb.bed
    j=$((j+3000000))
done
```

Next I’ll run bcftools view looping over each line in my partitions file


```python
#!/bin/bash

# Define the input BED file
bed_file="partitions_X_3mb.bed"

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

bed_file="partitions_X_3mb.bed"

module load python

# Loop through each VCF file matching the pattern

while IFS=$'\t' read -r chromosome start end; do
    # Define the output file prefix
    output_prefix="${chromosome}_${start}_${end}"

    for vcf_file in "${output_prefix}.vcf"; do
    python vcf2phylip.py -i "$vcf_file"
    done
done < "$bed_file"
```

Make trees for each alignment


```python
#!/bin/bash

bed_file="partitions_X_3mb.bed"

mkdir 3mbtrees_X

module load iqtree

# Loop through each VCF file matching the pattern

while IFS=$'\t' read -r chromosome start end; do
    # Define the output file prefix
    output_prefix="${chromosome}_${start}_${end}"

    # Loop through each VCF file matching the pattern
    for phy_file in "${output_prefix}.min4.phy"; do
    iqtree2 -s "$phy_file" -bb 1000 -wbt -nt AUTO 
    done

mv "$phy_file"* 3mbtrees_X

done < "$bed_file"
```

## Trees in sliding windows of 4MB across the X (partitioning vcf’s instead of the phylip)

Sliding window trees across X

Partitioning my vcf files into multiple files

First I’m going to make a bed file with information about each partition


```python
j=1
for i in {1..17}
do
    echo -e "ChrX_MullerAD\t${j}\t$((j+3999999))" >> partitions_X_4mb.bed
    j=$((j+4000000))
done
```

Next I’ll run bcftools view looping over each line in my partitions file


```python
#!/bin/bash

# Define the input BED file
bed_file="partitions_X_4mb.bed"

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

bed_file="partitions_X_4mb.bed"

module load python

# Loop through each VCF file matching the pattern

while IFS=$'\t' read -r chromosome start end; do
    # Define the output file prefix
    output_prefix="${chromosome}_${start}_${end}"

    for vcf_file in "${output_prefix}.vcf"; do
    python vcf2phylip.py -i "$vcf_file"
    done
done < "$bed_file"
```

Make trees for each alignment


```python
#!/bin/bash

bed_file="partitions_X_4mb.bed"

mkdir 4mbtrees_X

module load iqtree

# Loop through each VCF file matching the pattern

while IFS=$'\t' read -r chromosome start end; do
    # Define the output file prefix
    output_prefix="${chromosome}_${start}_${end}"

    # Loop through each VCF file matching the pattern
    for phy_file in "${output_prefix}.min4.phy"; do
    iqtree2 -s "$phy_file" -bb 1000 -wbt -nt AUTO 
    done

mv "$phy_file"* 4mbtrees_X

done < "$bed_file"
```
