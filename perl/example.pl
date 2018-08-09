#! /usr/bin/env perl
use Getopt::Long;
use Cwd;
use Time::Local;
use strict;
use warnings;

my $projecthome = "/home/rundata";
my $destination = $projecthome . "/";
my $top = getcwd();

my %projects = ( );

$projects{eeg_simJC}  = 1;
$projects{eeg_simstableJC}  = 1;



my ($help, $source, $scenario, $project);

GetOptions
(
		'source=s' => \$source, #ENVVARS assumed and added to say run22744_8_5_2008
        	'help' => \$help,
		'project=s' => \$project,
		
);

if(defined $help) {
	help();
	exit(0);
}

if(!(defined $project)) {
	$project = "eeg_simJC";
}

if(!(defined $source)) {
	print "Need directory to process!\n";
	die "Must use --source=xxxxxxxxx\n";
} else {
	if(!(-d "$source")) {
		print "Current files and folders here:\n";
		system("ls -l");
		die "Target dirrectory to process does not exist<$source>\n";
	}
}


#if(!(defined $scenario)) {
	#print "Need SCENARIO so we can name the new dataset\n";
	#die "Must use --scenario=xxxxxxxxx\n";
#}

#print "Pattern to find scenarios is <$scenario>\n";
#my @scenarios = `ls $scenario*`;


if(!(defined $project)) {
	print "Need PROJECT so we can place the new dataset\n";
	die "Must use --project=xxxxxxxxx\n";
} else {
	if(!(exists $projects{$project})) {
		print "Project choses:\n";
		foreach my $choice (keys %projects) {
			print "Project Choice: $choice\n";
		}
		die "Chosen project not in the list:<$project>\n";
	}
}


my $cwd = getcwd();

chdir($source);
my $source_wd = getcwd();


opendir SO, "." || die "Failed to pen dir  . ($source): !$\n";
	foreach my $file (readdir SO) {
		my  $thisdest = "";
		next if $file =~ /^\.\.?$/;
		if(!(-d $file)) {
			print "Skip regular file $file\n";
		}
		print "Consider $file\n";
		$thisdest = $destination . $project . "/" . "dataset_$file";


		if(-d $thisdest) {
			die "Dataset <$thisdest> already exists, chose new folder name <$file>\n";
		}

		system("mkdir $thisdest");
		if(!(-d $thisdest)) {
			die "Failed to create new dataset<$thisdest>\n";
		}

		system("mkdir $thisdest/shared");
		if(!(-d "$thisdest/shared")) {
			die "Failed to create new dataset<$thisdest/shared>\n";
		}

		#scope out the jobs in this set
		chdir("$file");

		my $dir = ".";

		opendir DH, $dir || die "Failed to pen dir $dir: !$\n";
		foreach my $targ (readdir DH) {
		
			next if $targ =~ /^\.\.?$/;
			print "looking at <$targ>\n";
			if($targ =~ /^.*num_subprocess.*$/) {
				print "Moving parallelism count targ<$targ>\n";
				print "cp $targ $thisdest/shared\n";
				system("cp $targ $thisdest/shared");
			 } elsif($targ =~ /^data.*?(\d+).mat$/) {
				mkdir("$thisdest/$1");
				print "cp $targ $thisdest/$1/input.mat\n";
				system("cp $targ $thisdest/$1/input.mat");
			 } elsif($targ =~ /^para.*\.mat$/) {
				print "Moving parameter targ<$targ>\n";
				system("cp $targ $thisdest/shared/param_file.mat");
			 } else {
			 	print "Ignoring <$targ>\n";
			 }
			
		
		}
		chdir("$source_wd");
		closedir(DH);
	}



closedir(SO);

print "Bye\n";
exit(0);

sub help {

		print "Usage:

./createdataset.pl --source=datain 	(folder with lors of jobs)



Options:
                [-h/--help]          See this
                [-p/--project=s]     (optional)What project is it for
                [-t/--source]        What directory to process
		 \n";

}
