grep "run complete" *.err | awk -F: '{print $1}'|sort > a
ls ../bam/*bam | grep -v temp | awk -F[/_] '{print $3".err"}' | sort > b 
diff <(grep "run complete" *.err | awk -F: '{print $1}'|sort) <(ls ../bam/*bam | grep -v temp | awk -F[/_] '{print $3".err"}' )
paste a b > c

