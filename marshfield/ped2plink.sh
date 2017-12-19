
## create a ped file in linkage format from the PLINK format
colrm 11 12 < test.ped > test2.ped

cut -d ' ' -f 6 test.ped > trait

paste test2.ped trait > test3.ped

## create data file in linkage format from the PLINK format

cut -d ' ' -f 2 test.map > test2.map

## open the test2.map and add one more line as trait

## create a ped file in linkage format from the PLINK format
colrm 11 12 < test.ped > test2.ped

cut -d ' ' -f 6 test.ped > trait

paste test2.ped trait > test3.ped

## create data file in linkage format from the PLINK format

cut -d ' ' -f 2 test.map > test2.map

## open the test2.map and add one more line as trait

paste pre_data test2.map > test3.map

## another way to edit the map file is to use R to read in the map file and only select the column with rs numbers as the second column, the first column is added as M as marker, and T as trait, the last line would be like 'T trait' save the file as test_link.dat

cp test3.map test_link.dat
cp test3.ped test_link.ped

## create bgl file from the above two files

java -jar /home/yez/WGI_Breast_Cancer/BEAGLE/BEAGLE_Utilities/linkage2beagle.jar test_link.dat test_link.ped > test_link.bgl

