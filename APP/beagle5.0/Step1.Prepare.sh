cd ~/hpc/db/hg19/beagle
for i in {1..22} X Y
do
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr$i.1kg.phase3.v5a.vcf.gz
done
wget http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/sample_info/20140625_related_individuals.txt
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/sample_info/integrated_call_male_samples_v3.20130502.ALL.panel
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/sample_info/integrated_call_samples.20130502.ALL.ped
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/sample_info/integrated_call_samples_v3.20130502.ALL.panel
unzip plink.GRCh37.map.zip
grep "EUR" integrated_call_male_samples_v3.20130502.ALL.panel | awk '{print $1}' > EUR.SampleList.txt
