```
cd /home/guosa/hpc/project/pmrp/phase1/plink
mkdir PA
plink --bfile FinalRelease_QC_20140311_Team1_Marshfield --allow-no-sex --pheno FinalRelease_QC_20140311_Team1_Marshfield.phen --pheno-name PheTyp3_PA_C1 --assoc counts --ci 0.95 --out ./PA/MCRI_PA_C2_Exom1
cp /home/guosa/hpc/project/pmrp/phase1/plink/PA/MCRI_PA_C2_Exom1.assoc 
cd /home/guosa/hpc/project/pmrp/phase2/PA
plink --bfile PMRP.PhaseII.Steven.Guo.PA.CEU --allow-no-sex --assoc counts --ci 0.95 --out MCRI_PA_C_Exom2
plink --meta-analysis /home/guosa/hpc/project/pmrp/phase1/plink/PA/MCRI_PA_C2_Exom1.assoc  /home/guosa/hpc/project/pmrp/phase2/PA/MCRI_PA_C1_Exom2.assoc  --out MCRI.PheTyp1_PA_C2_Exom1_Exom2
```
