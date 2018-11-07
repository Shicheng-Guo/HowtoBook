## Step 1:  Download Annotation Files
cd ~/hpc/tools/annovar
annotate_variation.pl -downdb -webfrom annovar -build hg19 dbnsfp35c  humandb/ &
cd -

annotate_variation.pl -downdb -webfrom annovar -build hg19 dbnsfp33a  humandb/
annotate_variation.pl -downdb -webfrom annovar -build hg19 ljb23_all  humandb/
annotate_variation.pl -downdb -webfrom annovar ljb23_fathmm humandb -buildver hg19  &
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar eigen humandb/
annotate_variation.pl -downdb -webfrom annovar -downdb avdblist humandb/ -build hg19 
annotate_variation.pl -downdb -webfrom annovar -build hg19 seq humandb/hg19_seq/
annotate_variation.pl -downdb wgEncodeGencodeBasicV19 humandb/ -build hg19
retrieve_seq_from_fasta.pl -format genericGene -seqdir ./humandb/hg19_seq/ -outfile ./humandb/hg19_wgEncodeGencodeBasicV19Mrna.fa ./humandb/hg19_wgEncodeGencodeBasicV19.txt 
annotate_variation.pl -downdb -webfrom annovar ljb23_metalr humandb/ -build hg19 &
annotate_variation.pl -downdb -webfrom annovar ljb23_gerp++ humandb/ -build hg19 &
annotate_variation.pl -downdb -webfrom annovar ljb23_all humandb/ -build hg19 &
annotate_variation.pl -downdb -webfrom annovar ljb23_metasvm humandb/ -build hg19 &
annotate_variation.pl -downdb wgEncodeGencodeBasicV19 humandb/ -build hg19
annotate_variation.pl -buildver hg19 -downdb wgEncodeGencodeBasicV19 humandb/
annotate_variation.pl -buildver hg19 -downdb dgvMerged humandb/
annotate_variation.pl -buildver hg19 -downdb tfbsConsSites humandb/
annotate_variation.pl -buildver hg19 -downdb wgRna humandb/
annotate_variation.pl -buildver hg19 -downdb targetScanS humandb/
annotate_variation.pl -buildver hg19 -downdb gwasCatalog humandb/
annotate_variation.pl -buildver hg19 -downdb cytoBand humandb/
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp147 humandb/
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp150 humandb/
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/ &
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar regsnpintron humandb/ &

# just for allele frequency
annotate_variation.pl -downdb -webfrom annovar exac03 humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar esp6500siv2 humandb -buildver hg19 &
annotate_variation.pl -downdb -webfrom annovar esp6500siv2_all humandb -buildver hg19 &
annotate_variation.pl -downdb -webfrom annovar gnomad_exome humandb -buildver hg19 &
# whole-exome data
annotate_variation.pl -downdb -webfrom annovar 1000g2015aug humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar kaviar_20150923 humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar hrcr1 humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar cg69 humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar gnomad_genome humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar dbnsfp30a humandb -buildver hg19 &
annotate_variation.pl -downdb -webfrom annovar esp6500siv2 humandb -buildver hg19 &
annotate_variation.pl -downdb esp6500siv2 humandb -buildver hg19 &
#  whole-genome data
cd ~/hpc/tools/annovar
annotate_variation.pl -downdb -webfrom annovar gerp++ humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar cadd humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar cadd13 humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar fathmm humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar eigen humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar gwava humandb -buildver hg19  &
# for CNV
annotate_variation.pl -downdb -webfrom annovar dbscsnv11 humandb -buildver hg19  
annotate_variation.pl -downdb -webfrom annovar spidex humandb -buildver hg19  
# disease-specific variants
annotate_variation.pl -downdb -webfrom annovar clinvar_20160302 humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar cosmic70 humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar icgc21 humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar nci60 humandb -buildver hg19  &
