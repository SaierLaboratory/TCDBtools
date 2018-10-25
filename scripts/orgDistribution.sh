#
#  Input file must be in NCBI standard FASTA format, where
#  headers start with a valid RefSeq accession. See the
#  following examples:
#
#       >CBN77049.1 PGPS/D6 [Ectocarpus siliculosus]
#       >CBN77049.1
#       >CBN77049
#
#  NOTE:
#  It is assumed that the NCBI NR and taxdb databases are locally
#  accessible via the blastdbcmd.
#


#Cleaning fasta headers
echo "cleaning fasta headers"
perl -i.bkp -pe 's/^\>(\w+).*$/\>$1/;' $1
echo

#Extracting accessions
echo "Extracting accessions"
grep '>' $1 | perl -pe 's/\>//;' > acc.txt
echo

#Extracting Superkingdom data
echo "Extracting Superkingdom data"
blastdbcmd -db nr -outfmt '%a %K' -entry_batch acc.txt -target_only -out kingdom.txt
echo

#Reporting superkingdom content
echo "Reporting superkingdom content"
cut -d " " -f 2 kingdom.txt | sort | uniq -c
