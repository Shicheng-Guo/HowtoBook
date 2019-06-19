#!/bin/bash
myvcf="$1"
# Print header
bcftools view -h -O v "$myvcf"
# Sort and remove duplicates per chromosome
for i in {1..22} MT X Y 
do
  bcftools view -H -r "$i" -O v "$myvcf" | sort -u -k2,2n -k4,4d -k5,5d
done
