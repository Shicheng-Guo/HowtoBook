## NIH Roadmap Epigenomics Project Data

### Data Download
Source: [WEB](https://www.ncbi.nlm.nih.gov/geo/roadmap/epigenomics/?view=matrix) [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE29611) [Web2](https://www.encodeproject.org/search/?type=Experiment&assay_title=WGBS&status=released&assembly=GRCh38&biosample_type=tissue&files.file_type=bigWig) [GEO](https://www.ncbi.nlm.nih.gov/geo/browse/?view=samples&display=200&series=16256&search=bisulfite%20sequencing&zsort=date)

Web2 is the best download source.

FTP: ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/roadmapepigenomics/by_experiment/

Wget: 
```
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE29nnn/GSE29611/suppl/GSE29611_RAW.tar
```

```md5sum
93ac8b8a17e1d561fdee5978a32b6dbe  BI.Adult_Liver.Bisulfite-Seq.3.hg19.bw
409d68942751b4bfe423910bd10bf09e  BI.Brain_Hippocampus_Middle.Bisulfite-Seq.149.hg19.bw
909f06da6978faa4c7f659e39542c250  BI.Brain_Hippocampus_Middle.Bisulfite-Seq.150.hg19.bw
8f4e4d962b4c4546df90fc8d6f609264  BI.Fetal_Muscle_Leg.Bisulfite-Seq.UW_H24996.hg19.bw
285276846af337eae652768010e2867b  BI.Fetal_Thymus.Bisulfite-Seq.UW_H24943.hg19.bw
6568ed396adfed031499f52bfa7c89ba  BI.hESC_Derived_CD184+_Endoderm_Cultured_Cells.Bisulfite-Seq.WGBS_Lib_19.hg19.bw
b0de5a3243f67b5e2f0be807748181a4  BI.hESC_Derived_CD184+_Endoderm_Cultured_Cells.Bisulfite-Seq.WGBS_Lib_6.hg19.bw
21d4c92f6caa8919d66f934e3a5baeb8  BI.hESC_Derived_CD56+_Ectoderm_Cultured_Cells.Bisulfite-Seq.WGBS_Lib_20.hg19.bw
a1503555caa79b271bc6860742020034  BI.hESC_Derived_CD56+_Ectoderm_Cultured_Cells.Bisulfite-Seq.WGBS_Lib_21.hg19.bw
a9a3a2b503a621116403e107a508d086  BI.hESC_Derived_CD56+_Ectoderm_Cultured_Cells.Bisulfite-Seq.WGBS_Lib_22.hg19.bw
3627d5216eaa53c2679c1532f7f83f8e  BI.hESC_Derived_CD56+_Ectoderm_Cultured_Cells.Bisulfite-Seq.WGBS_Lib_23.hg19.bw
ff29e2b1bdd1408468f94e2d29c42ec2  BI.hESC_Derived_CD56+_Mesoderm_Cultured_Cells.Bisulfite-Seq.WGBS_Lib_24.hg19.bw
dd865f9efe3e13390af1a255c38f7b9d  BI.hESC_Derived_CD56+_Mesoderm_Cultured_Cells.Bisulfite-Seq.WGBS_Lib_25.hg19.bw
4fde99efe44276844b65f0dc01237bc7  BI.HUES64.Bisulfite-Seq.WGBS_Lib_26.hg19.bw
f20075ea4a993c9292b1d7186d5504fe  BI.HUES64.Bisulfite-Seq.WGBS_Lib_39.hg19.bw
73d31caf048ce739a9f79d36a7335705  BI.Mobilized_CD34_Primary_Cells.Bisulfite-Seq.RO_01549.hg19.bw
89c1d2f59c9463a4e2deaf1be3668bcd  UCSD.Adipose_Tissue.Bisulfite-Seq.STL003.hg19.bw
a463290d3c0875b73a6864de084e007a  UCSD.Adrenal_Gland.Bisulfite-Seq.STL003.hg19.bw
2c4ec07cfca2fd9e775c7be02521625a  UCSD.Aorta.Bisulfite-Seq.STL003.hg19.bw
150764f4924f8cc3c54ff4a861a23475  UCSD.Esophagus.Bisulfite-Seq.STL003.hg19.bw
fd813d8907bd79310653f4aeb993f55c  UCSD.Gastric.Bisulfite-Seq.STL003.hg19.bw
a2ac6c3dca5a6f768a98aca3119e1a5f  UCSD.H1.Bisulfite-Seq.combined.hg19.bw
606e4c28ef73cd33b9fe737df467d394  UCSD.H1_BMP4_Derived_Mesendoderm_Cultured_Cells.Bisulfite-Seq.combined.hg19.bw
17bb51a15440fcbeec337f32e9bee24b  UCSD.H1_BMP4_Derived_Trophoblast_Cultured_Cells.Bisulfite-Seq.combined.hg19.bw
3fc6e0c75e1e6a5b6fafc041b48a0291  UCSD.H1_Derived_Mesenchymal_Stem_Cells.Bisulfite-Seq.combined.hg19.bw
af466113a99ba7ea7fc3e14bfe953875  UCSD.H1_Derived_Neuronal_Progenitor_Cultured_Cells.Bisulfite-Seq.combined.hg19.bw
48481ab9126ec3e1203c4bbaa852e3a5  UCSD.H9.Bisulfite-Seq.combined.hg19.bw
ab53fc1df2e064de9a38d961e59dd536  UCSD.IMR90.Bisulfite-Seq.combined.hg19.bw
7e928b399c2e963f73e7a411e677855a  UCSD.iPS_DF_19.11.Bisulfite-Seq.combined.hg19.bw
3910f20e05a9dee9a555ec62a48a9212  UCSD.iPS_DF_6.9.Bisulfite-Seq.combined.hg19.bw
276cc0e29a8440c825e189e6551df42a  UCSD.Left_Ventricle.Bisulfite-Seq.STL001.hg19.bw
33e8e115664655a671adcefc5d333785  UCSD.Left_Ventricle.Bisulfite-Seq.STL003.hg19.bw
34173b502b5e237d94b723b17f0843f4  UCSD.Lung.Bisulfite-Seq.STL002.hg19.bw
b985dca1c3521b44adcd56e3293cf49b  UCSD.Ovary.Bisulfite-Seq.STL002.hg19.bw
a112e3564c59cfa8658ce41410f041ba  UCSD.Pancreas.Bisulfite-Seq.STL003.hg19.bw
8d857e80333b0329f00d127d4be04fb0  UCSD.Psoas_Muscle.Bisulfite-Seq.STL003.hg19.bw
37e644f701bb5dd433b5564bad80775d  UCSD.Right_Atrium.Bisulfite-Seq.STL003.hg19.bw
0a11b6b8172a3e418bfb34c41090bbad  UCSD.Right_Ventricle.Bisulfite-Seq.STL003.hg19.bw
ef8edade291763a21d7b2de99e099520  UCSD.Sigmoid_Colon.Bisulfite-Seq.STL001.hg19.bw
68db3d30c028201dda89a704971036ee  UCSD.Sigmoid_Colon.Bisulfite-Seq.STL003.hg19.bw
16211648f712953d9e8aacfa66a36225  UCSD.Small_Intestine.Bisulfite-Seq.STL001.hg19.bw
79cffecf4797192bdd8278211a6367f7  UCSD.Spleen.Bisulfite-Seq.STL003.hg19.bw
3fcfd8303331ab095a93471ac400e4d5  UCSD.Thymus.Bisulfite-Seq.STL001.hg19.bw
626466c9840239a80611b1b4f2cfffb1  UCSF-UBC.Brain_Germinal_Matrix.Bisulfite-Seq.HuFGM02.hg19.bw
72cba19a5003fcc7b98758dbbf616547  UCSF-UBC.Breast_Luminal_Epithelial_Cells.Bisulfite-Seq.RM066.hg19.bw
50d6e2eaf85ee875246f4e14f10eec68  UCSF-UBC.Breast_Myoepithelial_Cells.Bisulfite-Seq.RM066.hg19.bw
caf52c836794834796079cd3566b76d3  UCSF-UBC.Neurosphere_Cultured_Cells_Cortex_Derived.Bisulfite-Seq.HuFNSC02.hg19.bw
bf14c7eb7fa62f37d46ba8de76fb907b  UCSF-UBC.Neurosphere_Cultured_Cells_Cortex_Derived.Bisulfite-Seq.HuFNSC04.hg19.bw
6a0c3bb20506de1932954acc8640058f  UCSF-UBC.Neurosphere_Cultured_Cells_Ganglionic_Eminence_Derived.Bisulfite-Seq.HuFNSC02.hg19.bw
7c47386d89dee8e1ad6920445c277e92  UCSF-UBC.Neurosphere_Cultured_Cells_Ganglionic_Eminence_Derived.Bisulfite-Seq.HuFNSC04.hg19.bw
cf80161415e4cd4a260dcc3c5655f078  UCSF-UBC.Penis_Foreskin_Fibroblast_Primary_Cells.Bisulfite-Seq.skin03.hg19.bw
4fbb20b73ec806c1d72532529947c9be  UCSF-UBC.Penis_Foreskin_Keratinocyte_Primary_Cells.Bisulfite-Seq.skin03.hg19.bw
49d2ec23be9985070ab38625bb5e1b1c  UCSF-UBC.Testis_Spermatozoa_Primary_Cells.Bisulfite-Seq.390ATA.hg19.bw
e5eaec6dda6e8588e3632fdd5df031d4  UCSF-UBC.UCSF-4star.Bisulfite-Seq.A21771.hg19.bw
5c9e58c56ad808463fc563eef0d72a5b  UCSF-UBC.UCSF-4star.Bisulfite-Seq.A21772.hg19.bw
072296cb2f10024038b69dd67a7a8ef6  UW.Fetal_Intestine_Large.Bisulfite-Seq.H-23769.lib-BiS.DS21177.hg19.bw
986073f488acbbe411468ee1feed814a  UW.Fetal_Intestine_Small.Bisulfite-Seq.H-23769.lib-BiS.DS21178.hg19.bw
```

### Data Enrollment
