#!/usr/bin/env bash
#
#Argument list:
#  1. RefSeq accession of the query protein
#  2. Accession of the protein in TCDB (Uniprot, trEMBL, RefSeq).
#  3. TCID of the system associated with accession passed to argument 2.
#
#Note: Accessions should not include the version number (e.g. .1, .2, etc.)
#

examineGBhit.pl -q $1 -s $2 -t $3 -o $1
