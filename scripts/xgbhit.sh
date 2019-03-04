#!/usr/bin/env bash
#
#Argument list:
#  1. RefSeq accession of the query protein
#  2. Accession of the protein in TCDB (Uniprot, trEMBL, RefSeq).
#  3. TCID of the system associated with accession passed to argument 2.
#  4. Substitution matrix to use (Optional. Default: BL62)
#Note: Accessions should not include the version number (e.g. .1, .2, etc.)
#


if [ "$1" = "-h" ]
then
  echo  
  echo "Arguments:"
  echo " 1  RefSeq accession of the query protein"
  echo ' 2. External Accession of the protein in TCDB (Uniprot, trEMBL, RefSeq).'
  echo ' 3. TCID of the system associated with accession passed to argument 2.'
  echo ' 4. Substitution matrix to use (Optional. Default: BL50)'
  exit 1
fi

#Define the substitution matrix to work with
mat="BL50"
if [ ! -z "$4" ]
then
  mat=$4
fi 
echo "Substitution Matrix: $mat"


examineGBhit.pl -q $1 -s $2 -t $3 -o $1 -m $mat
