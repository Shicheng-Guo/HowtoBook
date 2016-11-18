# how to get the time for the perl script with function of time

my $start_run = time();
my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job took $run_time seconds\n";

# how to get the time for the bash script with function of date
START=$(date +%s)
END=$(date +%s)
DIFF=$(echo "$END - $START" | bc)
echo "It takes DIFF=$DIFF seconds to complete this task..."
