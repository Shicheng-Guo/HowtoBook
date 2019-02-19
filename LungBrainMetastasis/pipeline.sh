/bin/BWA/bwa aln -f /output/FOO.sai \
            -t 3 /seq/REFERENCE/human_18.fasta \
            /seq/FQ/FOO.sanger.fq

/bin/BWA/bwa samse -f /output/FOO.sam \
            /seq/REFERENCE/human_18.fasta \
            /output/FOO.sai \
            /seq/FQ/FOO.sanger.fq

/bin/SAMTOOLS/samtools import /seq/REFERENCE/human_18.fasta \
            /output/FOO.sam \
            /output/FOO.bam

/bin/SAMTOOLS/samtools sort /output/FOO.bam /output/FOO.sorted

/bin/SAMTOOLS/samtools index /output/FOO.sorted.bam \
            /output/FOO.sorted.bam.bai

java -jar /bin/GTK/GenomeAnalysisTK.jar -T RealignerTargetCreator \
            -R /seq/REFERENCE/human_18.fasta \
            -I /output/FOO.sorted.bam \
            -o /output/FOO.intervals

java -jar /bin/GTK/GenomeAnalysisTK.jar -T IndelRealigner 
            -R /seq/REFERENCE/human_18.fasta \
            -I /output/FOO.sorted.bam \
            -targetIntervals /output/FOO.intervals \
            --output /output/FOO.sorted.realigned.bam

/bin/SAMTOOLS/samtools index /output/FOO.sorted.realigned.bam \
            /output/FOO.sorted.realigned.bam.bai
java -jar /bin/GTK/GenomeAnalysisTK.jar -T IndelGenotyperV2 \
            -R /seq/REFERENCE/human_18.fasta \
            -I /output/FOO.sorted.realigned.bam \
            -O /output/FOO_indel.txt \
            --verbose \
            -o /output/FOO_indel_statistics.txt

java -jar /bin/GTK/GenomeAnalysisTK.jar -T UnifiedGenotyper \
            -R /seq/REFERENCE/human_18.fasta \
            -I /output/FOO.sorted.realigned.bam \
            -varout /output/FOO.geli.calls \
            -vf GELI \
            -stand_call_conf 30.0 \
            -stand_emit_conf 10.0 \
            -pl SOLEXA
