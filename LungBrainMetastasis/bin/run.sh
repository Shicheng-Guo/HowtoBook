perl VEPmerge.pl > BLM-VEP.hg19.bed
perl bed2mutationprofile.pl > LBM.MutationProfile.txt
cp LBM.MutationProfile.txt LBM.MutationProfile.test.txt
perl -p -i -e 's/3_prime_UTR_variant/regulatory_region_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/5_prime_UTR_variant/regulatory_region_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/coding_sequence_variant/missense_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/downstream_gene_variant/Others/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/frameshift_variant/frameshift_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/inframe_deletion/frameshift_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/intergenic_variant/Others/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/intron_variant/Others/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/mature_miRNA_variant/regulatory_region_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/missense_variant/missense_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/NMD_transcript_variant/missense_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/non_coding_transcript_exon_variant/regulatory_region_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/non_coding_transcript_variant/regulatory_region_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/protein_altering_variant/missense_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/regulatory_region_variant/regulatory_region_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/splice_acceptor_variant/splice_related_variants/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/splice_donor_variant/splice_related_variants/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/splice_region_variant/splice_related_variants/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/start_lost/missense_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/stop_gained/missense_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/stop_lost/missense_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/synonymous_variant/Others/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/TF_binding_site_variant/regulatory_region_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/upstream_gene_variant/regulatory_region_variant/g' LBM.MutationProfile.test.txt
