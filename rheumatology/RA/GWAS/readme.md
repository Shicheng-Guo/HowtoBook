Genetic Association SNPs from GWAS-Meta analysis (20190208)

[128 SNPs](GWAS-Meta-128-SNPs.20190208.vcf) have been identified in previous meta-analysis to GWAS-RA data. However, we can find there are more than GWAS-Catalog RA-SNPs were recorded. 


```
Response to TNF antagonist treatment
Anti-cyclic Citrullinated Peptide Antibody
Response to tocilizumab in rheumatoid arthritis
Response to TNF-alpha inhibitors in rheumatoid arthritis
Response to anti-TNF therapy in rheumatoid arthritis
Response to methotrexate in rheumatoid arthritis
Rheumatoid arthritis (ACPA-negative)
Joint damage progression in ACPA-negative rheumatoid arthritis
Response to anti-TNF therapy in rheumatoid arthritis
Rheumatoid arthritis (rheumatoid factor and/or anti-cyclic citrullinated peptide seropositive)
Rheumatic heart disease
Bone erosion in rheumatoid arthritis
Severe progression in rheumatoid arthritis
Juvenile idiopathic arthritis (oligoarticular or rheumatoid factor-negative polyarticular)
Sero-negative rheumatoid arthritis
ACPA-positive rheumatoid arthritis (smoking interaction)
Response to TNF inhibitor in rheumatoid arthritis (change in swollen 28-joint count)
```

UK-biobank data re-analysis
```
awk -F'[:|\t]' '{print $1,$2,$3,$4,$5,$6,$9,$10,$12}' OFS="\t" 20002_1464.assoc.tsv  > 20002_1464.assoc.txt &
awk -F'[:|\t]' '{print $1,$2,$3,$4,$5,$6,$9,$10,$12}' OFS="\t" 40001_M069.assoc.tsv  > 40001_M069.assoc.txt  &
awk -F'[:|\t]' '{print $1,$2,$3,$4,$5,$6,$9,$10,$12}' OFS="\t" M05.assoc.tsv  > M05.assoc.txt &
awk -F'[:|\t]' '{print $1,$2,$3,$4,$5,$6,$9,$10,$12}' OFS="\t" M06.assoc.tsv  > M06.assoc.txt &

echo "CHRBP\tA1\tA2\tSNP\tN\tBETA\tSE\tP" > head.txt
cat head.txt 20002_1464.assoc.txt > 20002_1464.assoc.hg19.txt &
cat head.txt 40001_M069.assoc.txt > 40001_M069.assoc.hg19.txt &
cat head.txt M05.assoc.txt > M05.assoc.hg19.txt &
cat head.txt M06.assoc.txt > M06.assoc.hg19.txt &

grep -v 'variant' 20002_1464.assoc.hg19.txt > 20002_1464.assoc.hg19.trim.txt &
grep -v 'variant' 40001_M069.assoc.hg19.txt > 40001_M069.assoc.hg19.trim.txt &
grep -v 'variant' M05.assoc.hg19.txt > M05.assoc.hg19.trim.txt &
grep -v 'variant' M06.assoc.hg19.txt > M06.assoc.hg19.trim.txt &

plink --meta-analysis 40001_M069.assoc.hg19.trim.txt 20002_1464.assoc.hg19.trim.txt M05.assoc.hg19.trim.txt M06.assoc.hg19.trim.txt  + logscale  --out UK-biobank-RA
```
