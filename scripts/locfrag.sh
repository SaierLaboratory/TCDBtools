
#Help
if [ "$1" = "-h" ]
then
  echo "Locate aligned fragments within 2 full proteins"
  echo "Arguments:"
  echo " 1. NCBI Accession of first protein 1"
  echo " 2. Aligned fragment of protein 1" 
  echo " 3. NCBI Accession of first protein 2"
  echo " 4. Aligned fragment of protein 2"
  echo " 5. (Optional) substitution matrix to use. (Defaul: BL62)"
  exit 1
fi



#Define the substitution matrix to work with
mat="BL50"
if [ ! -z "$5" ]
then
  mat=$5
fi 
echo "Substitution Matrix: $mat"



locateFragment.pl -a  $1 -f $2 
locateFragment.pl -a  $3 -f $4
alignSeqsFiles.pl -q $1_frag.faa -ql $1_frag -s $3_frag.faa -sl $3_frag -e 10 -c 20 -cc X -m $mat
alignSeqsFiles.pl -q $1.faa -ql $1 -s $3.faa -sl $3 -e 10 -c 5 -cc X -m $mat
open ssearch*/*.html
