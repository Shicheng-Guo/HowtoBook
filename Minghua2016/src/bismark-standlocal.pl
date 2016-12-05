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
 
my $job_file_name = "bismark.sh";
open(OUT, ">$job_file_name") || die("Error in opening file $job_file_name.\n");   
    
open F,"$sample_group_file";
while(<F>){
    chomp;
    next if /^\s+$/;
    my @read = split /\t/;
    my ($sample,undef)=split /_1.fastq.gz/,$read[0];
    my $project="bismark";
    my $analysis="alignment";
    
    my $curr_dir = $dir;
    
    my $ppn=8;
    my $multicore=2;
    my $walltime="167:00:00";
    my $queue="hotel"; # hotel,pdafm,condo
    # chomp(my $phredcheck=`perl /home/shg047/bin/checkphred.pl $read[0]`);
    # my ($phred)=split /\s+/,$phredcheck;
    if(scalar(@read) eq 2){
    my($sample,undef)=split /_R1.fastq.gz/,$read[0]; 	
    my($sample1,undef)=split /.fastq.gz/,$read[0];
    my($sample2,undef)=split /.fastq.gz/,$read[1];
        print OUT "cd $dir\n";
	#print OUT "fastq-dump --split-files --gzip $SRR\n";
#	print OUT "trim_galore --paired --phred$phred --fastqc --illumina $sample1.fastq.gz $sample2.fastq.gz --output_dir ../fastq_trim\n";
#	print OUT "bismark --bowtie2 --multicore 4 --phred$phred-quals --fastq -L 21 -N 1 ~/db/hg19/align/bismark/ -1 ../fastq_trim/$sample1\_val_1.fq.gz -2 ../fastq_trim/$sample2\_val_2.fq.gz -o ../bam\n";
#	print OUT "filter_non_conversion --paired ../bam/$sample1\_val_1_bismark_bt2_pe.bam\n";
#	print OUT "samtools sort ../bam/$sample1\_val_1_bismark_bt2_pe.nonCG_filtered.bam ../sortbam/$sample\_bismark_bt2_pe.sort\n";
#	print OUT "samtools index ../sortbam/$sample\_bismark_bt2_pe.sort.bam\n";
	print OUT "bismark_methylation_extractor --no_overlap --merge_non_CpG --cutoff 10 --multicore 8 --paired-end --bedGraph --ignore 1 --buffer_size 4G --zero_based --comprehensive --output ../methyfreq  ../bam/$sample1\_val_1_bismark_bt2_pe.nonCG_filtered.bam\n";
        print OUT "bedtools intersect -wa -a ../methyfreq/$sample1\_val_1_bismark_bt2_pe.nonCG_filtered.bedGraph.gz.bismark.zero.cov -b ../bed/target.bed > ../bedgraph/$sample.bedgraph\n\n";
    }else{
    my($sample,undef)=split /.fastq.gz/,$read[0]; 	
    my($sample1,undef)=split /.fastq.gz/,$read[0];
	#print OUT "fastq-dump --split-files --gzip $SRR\n";
	print OUT "trim_galore --phred$phred --fastqc --illumina $sample1.fastq.gz --output_dir ../fastq_trim\n";
    print OUT "bismark --bowtie2 --non_directional $phred-quals --fastq -L 32 -N 1 --multicore 6 ~/db/hg19/align/bismark/ ../fastq_trim/$sample1\_trimmed.fq.gz -o ../bam\n";  
	print OUT "samtools sort ../bam/$sample\_se.bam -o ../sortbam/$sample.sort.bam\n";
	print OUT "samtools index ../sortbam/$sample.sort.bam\n";
	print OUT "bismark_methylation_extractor --mapq 30 --no_overlap --single-end --bedGraph --ignore 3 --buffer_size 4G --zero_based --comprehensive --output ../methyfreq  ../bam/$sample\_se.bam";
    }
}
close(OUT);
# filter_non_conversion --paired ../bam/CU1_R1_val_1_bismark_bt2_pe.bam
# samtools sort ../bam/CU1_R1_val_1_bismark_bt2_pe.nonCG_filtered.bam -o ../sortbam/CU1_bismark_bt2_pe.bam
# samtools index ../sortbam/CU1_bismark_bt2_pe.bam
# bismark_methylation_extractor --single-end --bedGraph --ignore 3 --buffer_size 4G --zero_based --comprehensive --output ../methyfreq  ../bam/CU1_R1_val_1_bismark_bt2_pe.nonCG_removed_seqs.bam


















