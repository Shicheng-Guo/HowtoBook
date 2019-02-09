cd /gpfs/home/guosa/hpc/db/hg19/1000Genome

mkdir temp
for i in {1..22} X Y
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo tabix -p vcf ALL.chr$i.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz >> $i.job
qsub $i.job
done




#######################################################
AHRR: chr5:317989-325025
cd /home/guosa/hpc/methylation/methtype
plink --bfile /home/guosa/hpc/db/hg19/1000Genome/chr5 --chr 5 --from-bp 317989 --to-bp 325025 --recode fastphase --freq --out /home/guosa/hpc/methylation/methtype/AHRR
plink --bfile /home/guosa/hpc/db/hg19/1000Genome/chr5 --chr 5 --from-bp 317989 --to-bp 325025 --recode fastphase --freq --out /home/guosa/hpc/methylation/methtype/AHRR

cd /home/guosa/hpc/db/hg19/plan2/commonCpGSNP
rm cpgSNP.hg19.uni.bed
for i in `ls *cpgsnp.bed`
do
perl unitrim.pl $i >> cpgSNP.hg19.uni.bed
done
awk '{print $4}' cpgSNP.hg19.uni.bed > cpgSNP.hg19.list.txt 

cd /home/guosa/hpc/cpgSNP
for i in {1..22} X Y
do
plink --bfile /home/guosa/hpc/db/hg19/1000Genome/chr$i --extract /home/guosa/hpc/db/hg19/plan2/commonCpGSNP/cpgSNP.hg19.list.txt --make-bed --out chr$i.cpgSNP
done

cd /home/guosa/hpc/cpgSNP
for i in {1..22} X Y
do
perl assignRefence.pl chr$i.cpgSNP.bim > chr$i.cpgSNP.bim.tr
mv chr$i.cpgSNP.bim.tr chr$i.cpgSNP.bim
awk '{print $2"\tC"}' chr$i.cpgSNP.bim > chr$i.C.ref
plink --bfile chr$i.cpgSNP --reference-allele cpgSNP.ref.txt --make-bed --out chr$i.cpgSNP.C
done

plink --bfile /home/guosa/hpc/cpgSNP/chr5.cpgSNP.C --chr 5 --from-bp 317989 --to-bp 325025 --recode fastphase --reference-allele cpgSNP.ref.txt  --freq --out /home/guosa/hpc/cpgSNP/AHRR
cd /home/guosa/hpc/cpgSNP
perl swith2csv.pl fastphase_hapguess_switch.out > AHRR.hap.txt

plink --bfile FGF6 --chr 12 --from-bp 4536137 --to-bp 4561951 --make-bed --freq counts --allow-no-sex --out data2
plink --bfile FGF6 --chr 12 --from-bp 4536137 --to-bp 4561951 --make-bed --freq case-control --allow-no-sex --out FGF6-C2
data<-read.table("AHRR.txt",row.names=1)
colSums(input)[colSums(input)>50]
awk '{print $1}' two_alof_all_combed_v4.phen > two.phen.txt
awk '{print $1}' FGF6.fam > one.phen.txt
diff one.phen.txt two.phen.txt

MC.154314@1075641920 MC.154314@1075641920 0 0 0 -9
MC.154316@1075680154 MC.154316@1075680154 0 0 0 -9

data<-read.table("two_alof_all_combed_v4.phen",head=T)
output<-data.frame(FID=data$FID,FID=data$FID,0,0,0,PheTyp7_Iron_C1=data$PheTyp7_Iron_C1)
write.table(output,file="FGF6-C1.phen",sep="\t",quote=F,col.names=F,row.names=F)
data<-read.table("FGF6.dip.txt",head=F,row.names=1)
saminfo<-read.table("two_alof_all_combed_v4.phen",head=T)

levels(data[,1])[4]=0
levels(data[,2])[3]=0
levels(data[,1])[2]="CTCGGCCCAC"
levels(data[,1])[7]="M1"
levels(data[,2])[6]="M1"
levels(data[,1])[5]="M4"
levels(data[,2])[4]="M4"

levels(data[,1]) %in% levels(data[,2])
levels(data[,1])[2]<-"CTCGGCCCAC"
table(data[,1],data[,2])


write.table(table(data[,1],data[,2]),file="FGF6.PMRP.txt",sep="\t",col.names=NA,row.names=T,quote=F)

CON1<-which(saminfo$PheTyp7_Iron_C1==1)
CAS1=which(saminfo$PheTyp7_Iron_C1==2)
write.table(table(data[CAS1,1],data[CAS1,2]),file="FGF6.PMRP.C1.Case.tsv",sep="\t",col.names=NA,row.names=T,quote=F)
write.table(table(data[CON1,1],data[CON1,2]),file="FGF6.PMRP.C1.CONTROL.tsv",sep="\t",col.names=NA,row.names=T,quote=F)

CON2=which(saminfo$PheTyp7_Iron_C2==1)
CAS2=which(saminfo$PheTyp7_Iron_C2==2)
write.table(table(data[CAS2,1],data[CAS2,2]),file="FGF6.PMRP.C2.Case.txt",sep="\t",col.names=NA,row.names=T,quote=F)
write.table(table(data[CON2,1],data[CON2,2]),file="FGF6.PMRP.C2.CONTROL.txt",sep="\t",col.names=NA,row.names=T,quote=F)

tab1<-table(data[CAS2,1],data[CAS2,2])
tab2<-table(data[CON2,1],data[CON2,2])

table(C3)
table(C4)
C4[C4$V2=="M1",]


/home/guosa/hpc/wgbs/GSE17972
/home/guosa/hpc/methylation/methtype

for i in `ls *.cov.gz`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -q longq  >> $i.job
echo cd $(pwd) >> ${i}.job
echo gunzip $i >> ${i}.job
qsub ${i}.job
done


mkdir ../mf
for i in `ls *.cov.gz`
do
j=${i/_L001_R1_001_00_bismark_bt2_pe.bam/}
BismarkRefereDb="~/hpc/db/hg19/bismrk/"
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=32 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -q longq  >> $i.job
echo cd $(pwd) >> ${i}.job
echo samtools sort -n $i $i.nsort >> ${i}.job
echo $j
done


wigToBigWig GSM916049_BI.Adult_Liver.Bisulfite-Seq.3.wig ~/hpc/db/hg19/hg19.chrom.sizes GSM916049_BI.Adult_Liver.Bisulfite-Seq.3.bw 
bigWigAverageOverBed GSM916049_BI.Adult_Liver.Bisulfite-Seq.3.bw ~/hpc/db/hg19/BUR.GRCH37.bed GSM916049_BI.Adult_Liver.Bisulfite-Seq.3.bw.tab




