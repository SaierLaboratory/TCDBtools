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
  echo " 1  RefSeq accession of the query protein (mandatory)"
  echo ' 2. External Accession of the protein in TCDB (i.e. Uniprot, trEMBL, RefSeq) (mandatory)'
  echo ' 3. TCID of the system associated with accession passed to argument 2 (mandatory).'
  echo ' 4. BlastDB to extract the sequence from the query protein (optional).'
  echo ' 5. Substitution matrix to use (Optional. Default: BL50)'
  exit 1
fi

#Define the substitution matrix to work with
mat="BL50"
bdb='nr'
if [[ ! -z "$4" ]] && ([[ "$4" == "BL50" ]] || [[ "$4" == "BL62" ]])
then
    mat=$4
elif [[ ! -z "$4" ]]
then
    bdb=$4    
fi

if [[ ! -z "$5" ]]
then
    mat=$5
fi


echo "Substitution Matrix: $mat"
echo "BlastDB: $bdb"


examineGBhit.pl -q $1 -s $2 -t $3 -o $1 -m $mat -bdb $bdb
