#!/usr/bin/perl -w

# bismark to alignment single and pair-end fastq in same project
# Good habit to creat parameter table for a project
# Contact: Shicheng Guo
# Version 1.3
# Update: Nov/23/2016

use strict;
use Cwd;
my $dir=getcwd;
die USAGE() if scalar(@ARGV<2);

my $sample_group_file=shift @ARGV;
my $phred=shift @ARGV;
my $submit=shift @ARGV;

my $humanBismarkRefereDb="/home/shg047/db/hg19/bismark/";
my $mouseBismarkRefereDb="/home/shg047/db/mm9/bismark/";

mkdir "../fastq_trim" if ! -e "../fastq_trim";
mkdir "../bam" if ! -e "../bam";
mkdir "../bedgraph" if ! -e "../bedgraph";
mkdir "../sortbam" if ! -e "../sortbam";
mkdir "../methyfreq" if ! -e "../methyfreq";
 

open F,"$sample_group_file";
while(<F>){
    chomp;
    next if /^\s+$/;
    my @read = split /\t/;
    my ($sample,undef)=split /_1.fastq.gz/,$read[0];
    my $project="bismark";
    my $analysis="alignment";
    my $job_file_name = $sample . ".job";
    my $curr_dir = $dir;
    
    my $ppn=1;
    my $multicore=1;
    my $walltime="8:00:00";
    my $queue="condo"; # hotel,pdafm,condo
    # chomp(my $phredcheck=`perl /home/shg047/bin/checkphred.pl $read[0]`);
    # my ($phred)=split /\s+/,$phredcheck;
    
    open(OUT, ">$job_file_name") || die("Error in opening file $job_file_name.\n");   
    print OUT "#!/bin/csh\n";
    print OUT "#PBS -N $job_file_name\n";
    print OUT "#PBS -q $queue\n";  # glean is free
    print OUT "#PBS -l nodes=1:ppn=$ppn\n";
    print OUT "#PBS -l walltime=$walltime\n";
    print OUT "#PBS -o $sample.log\n";
    print OUT "#PBS -e $sample.err\n";
    print OUT "#PBS -V\n";
    print OUT "#PBS -M shicheng.guo\@gmail.com \n";
    print OUT "#PBS -m abe\n";
    print OUT "#PBS -A k4zhang-group\n";
    print OUT "cd $curr_dir\n";
    print "$job_file_name\n";
    if(scalar(@read) eq 2){
    my($sample,undef)=split /_R1.fastq.gz/,$read[0]; 	
    my($sample1,undef)=split /.fastq.gz/,$read[0];
    my($sample2,undef)=split /.fastq.gz/,$read[1];
	#print OUT "fastq-dump --split-files --gzip $SRR\n";
	print OUT "trim_galore --paired --phred$phred --fastqc --illumina $sample1.fastq.gz $sample2.fastq.gz --output_dir ../fastq_trim\n";
	print OUT "bismark --bowtie2 --non_directional --phred$phred-quals --fastq -L 21 -N 1 /home/shg047/db/hg19/bismark/ -1 ../fastq_trim/$sample1\_val_1.fq.gz -2 ../fastq_trim/$sample2\_val_2.fq.gz -o ../bam\n";
	print OUT "filter_non_conversion --paired ../bam/$sample1\_val_1_bismark_bt2_pe.bam\n";
	print OUT "samtools sort ../bam/$sample1\_val_1_bismark_bt2_pe.nonCG_filtered.bam -o ../sortbam/$sample\_bismark_bt2_pe.sort.bam\n";
	print OUT "samtools index ../sortbam/$sample\_bismark_bt2_pe.sort.bam\n";
	print OUT "bismark_methylation_extractor --single-end --bedGraph --ignore 3 --buffer_size 4G --zero_based --comprehensive --output ../methyfreq  ../bam/$sample1\_val_1_bismark_bt2_pe.nonCG_filtered.bam";
    }else{
    my($sample,undef)=split /.fastq.gz/,$read[0]; 	
    my($sample1,undef)=split /.fastq.gz/,$read[0];
	#print OUT "fastq-dump --split-files --gzip $SRR\n";
	print OUT "trim_galore --phred$phred --fastqc --illumina $sample1.fastq.gz --output_dir ../fastq_trim\n";
    print OUT "bismark --bowtie2 --non_directional $phred-quals --fastq -L 32 -N 1 --multicore 6 /home/shg047/db/hg19/bismark/ ../fastq_trim/$sample1\_trimmed.fq.gz -o ../bam\n";  
	print OUT "samtools sort ../bam/$sample\_se.bam -o ../sortbam/$sample.sort.bam\n";
	print OUT "samtools index ../sortbam/$sample.sort.bam\n";
	print OUT "bismark_methylation_extractor --single-end --bedGraph --ignore 3 --buffer_size 4G --zero_based --comprehensive --output ../methyfreq  ../bam/$sample\_se.bam";
    }
   close(OUT);
   if($submit eq 'submit'){
   system("qsub $job_file_name");
   }
}

# filter_non_conversion --paired ../bam/CU1_R1_val_1_bismark_bt2_pe.bam
# samtools sort ../bam/CU1_R1_val_1_bismark_bt2_pe.nonCG_filtered.bam -o ../sortbam/CU1_bismark_bt2_pe.bam
# samtools index ../sortbam/CU1_bismark_bt2_pe.bam
# bismark_methylation_extractor --single-end --bedGraph --ignore 3 --buffer_size 4G --zero_based --comprehensive --output ../methyfreq  ../bam/CU1_R1_val_1_bismark_bt2_pe.nonCG_removed_seqs.bam



sub USAGE{
print "\nUSAGE: $0 Samconfigfile phredscore submit\n";
print '
--------------------------------------------------------------------------------------------------------------
This command will trim the Fastq files with trim_galore for single-end and pair-end BS-seq data simultaniously. 
Example:  perl fastq2trimfastq.pl FastqMatchConfig.txt submit
SamConfigFile Format:
HOT243.fq.gz
HOT265.fq.gz
BMT1.read1.fq.gz        BMT1.read2.fq.gz
BMT2.read1.fq.gz        BMT2.read2.fq.gz
submit|notsumbit: indicate qsub job to TSCC or not?
The trimed fastq will be save in ../fastq_trim. 
Previous Script: 
1, Download or Prepare SRA Download Configure  (http://www.ebi.ac.uk/ena/data/view/SRP028600)
2, perl fastqDownload.pl SamConfig.txt 8 submit   (fastq-dump --split-files --gzip)
3, trim_galore --phred33 --fastqc --stringency 3 -q 20 --trim1 --length 20 --gzip --clip_R1 2 --three_prime_clip_R1 2 --illumina SRR299055_1.fastq.gz --output_dir ../fastq_trim
4, bismark_genome_preparation ./     # merge all the fa to one file (mm9.fa or hg19.fa)
--------------------------------------------------------------------------------------------------------------
';


}















