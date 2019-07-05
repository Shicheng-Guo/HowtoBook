vcf="$1"
bcftools view -h -O v "$vcf"
for i in {1..22} MT X Y 
do
bcftools view -H -r "$i" -O v "$vcf" | sort -u -k2,2n -k4,4d -k5,5d
done

# Run like this way:
# sh remove_VCF_duplicates.sh All_samples_Exome.vcf.gz > All_samples.undup.vcf
