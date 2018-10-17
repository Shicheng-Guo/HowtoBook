#!/bin/bash
if [ ! -f beagle.28Sep18.793.jar ]; then
  echo
  echo "Downloading beagle.28Sep18.793.jar"
  wget http://faculty.washington.edu/browning/beagle/beagle.28Sep18.793.jar
fi

if [ ! -f bref3.28Sep18.793.jar ]; then
  echo
  echo "Downloading bref3.28Sep18.793.jar"
  wget http://faculty.washington.edu/browning/beagle/bref3.28Sep18.793.jar
fi

echo

if [ ! -f test.28Sep18.793.vcf.gz ]; then
    echo
    echo "*** Downloading some 1000 Genomes Project data to file: test.28Sep18.793.vcf.gz ***"
    wget http://faculty.washington.edu/browning/beagle/test.28Sep18.793.vcf.gz
fi

echo
echo "*** Creating test files: ref.28Sep18.793.vcf.gz target.28Sep18.793.vcf.gz ***"
echo
zcat test.28Sep18.793.vcf.gz | cut -f1-190 | tr '/' '|' | gzip > ref.28Sep18.793.vcf.gz
zcat test.28Sep18.793.vcf.gz | cut -f1-9,191-200 | gzip > target.28Sep18.793.vcf.gz

echo
echo "*** Running test analysis with \"gt=\" argument ***"
echo
java -jar beagle.28Sep18.793.jar gt=test.28Sep18.793.vcf.gz out=out.gt

echo
echo "*** Running test analysis with \"ref=\" and \"gt=\" arguments ***"
echo
java -jar beagle.28Sep18.793.jar ref=ref.28Sep18.793.vcf.gz gt=target.28Sep18.793.vcf.gz out=out.ref

echo
echo "*** Making \"bref3\" file ***"
echo
java -jar bref3.28Sep18.793.jar ref.28Sep18.793.vcf.gz > ref.28Sep18.793.bref3

echo
echo "*** Running test analysis with \"bref3\" file ***"
echo
java -jar beagle.28Sep18.793.jar ref=ref.28Sep18.793.bref3 gt=target.28Sep18.793.vcf.gz out=out.bref3
