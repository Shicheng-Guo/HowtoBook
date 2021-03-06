## Genome-wide DNA methylation analysis to Roadmap and Encode project

## Data
Download WGBS dataset: [42 WGBS dataset](https://www.encodeproject.org/search/?type=Experiment&status=released&assay_slims=DNA+methylation&biosample_type=tissue&y.limit=&assay_title=WGBS&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&replicates.library.biosample.life_stage=fetal&replicates.library.biosample.life_stage=adult&replicates.library.biosample.life_stage=child&files.file_type=fastq&files.file_type=bed+bedMethyl&files.file_type=bigBed+bedMethyl&files.file_type=bam&files.file_type=bigWig&files.file_type=sra&files.run_type=single-ended&lab.title=Bradley+Bernstein%2C+Broad&lab.title=Joe+Ecker%2C+Salk&replication_type=unreplicated)
```bash 
xargs -n 1 curl -O -L < files.txt
```

```bash
wget -r ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/roadmapepigenomics/by_experiment/Bisulfite-Seq
awk '{print $1,$2,$3,$1":"$2"-"$3}' OFS="\t" wgEncodeRegTfbsClusteredV3.bed > wgEncodeRegTfbsClustered.bed
cd /home/shg047/oasis/Roadmap/wig
```
