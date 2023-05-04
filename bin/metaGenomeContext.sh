#!/usr/bin/env bash

#Input file must be the  the file in the assembly ending with:
#  _feature_table.txt.gz

if [ "$1" = "-h" ]
then
  echo "Find the genomic context of a gene within a meta genome"
  echo "Arguments:"
  echo " 1. File from assembly ending with: _feature_table.txt.gz"
  echo " 2. Unique accession of a gene in genome (e.g. RefSeq or Locus Tag" 
  exit 1
fi


AssFile=${1?Assembly ID is mandatory}
Acc=${2?Protein accession is mandatory}



function mgc () {

  #First find the Assembly piece to which the gene belongs
  FragmentID=$(zcat $1 | grep $2 | cut -f 7) 

  #now get all the genes in that assembly piece
  zgrep $FragmentID $1 | grep CDS | cut -f 7-11,14 | sort -k2n

}

mgc $AssFile $Acc

