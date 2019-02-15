Whole-exome sequencing identifies somatic mutations associated with mortality and lung cancer metastasis to the brain



Mathod and Materials

* 100ng input DNA (cancer and normal), 200-350bp fragment, SureSelect v2 Exome bait (Agilent), HiSeq, average exome coverage of 83.3x 
* Phylogenetic Tree Construction: All nonsilent mutations that passed validation (EAC001, EAC003,and EAC005) or further fi ltering (EAC006, EAC009, EAC014, EAC015,and EAC017) were considered for the purpose of determining phylogenetic trees. Trees were built using binary presence/absence matrices built from the regional distribution of variants within the tumor. The R Bioconductor package phangorn was utilized to perform the parsimony ratchet method, generating unrooted trees. Branch lengths were determined using the acctran function
* Canopy, a method for inferring the evolutionary phylogeny of a tumor using both somatic copy number alterations and single-nucleotide
alterations from one or more samples derived from a single patient. Canopy is applied to bulk sequencing datasets of both longitudinal and spatial experimental designs and to a transplantable metastasis model derived from human cancer cell line MDA-MB-231

