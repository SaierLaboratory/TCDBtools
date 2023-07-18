
#Help
if [ "$1" = "-h" ]
then
  echo "Locate aligned fragments within 2 full proteins"
  echo "Arguments:"
  echo " 1. Accession of protein 1"
  echo " 2. Aligned fragment of protein 1" 
  echo " 3. Accession of protein 2"
  echo " 4. Aligned fragment of protein 2"
  echo " 5. (optional) blast DB to extract sequences from (Default: uniref90)"
  echo " 6. (Optional) substitution matrix to use. (Defaul: BL50)"
  echo " 7. (Optional) Indicate whether plots will be shown (Values: show/quiet; Default: show)" 
  exit 1
fi



#Define the substitution matrix to work with
mat="BL50"
mode=""
db="uniref90"

#Identify the blast DB to extract sequences from
if [[ ! -z "$5" ]] && ([[ "$5" != "quiet" ]] && [[ "$5" != "show" ]] && [[ "$5" != "uniref90" ]])
then
  db="$5"
fi 

#Identify the type of substitution matrix, if given
if [[ ! -z "$6" ]] && ([[ "$6" != "quiet" ]] && [[ "$6" != "show" ]] && [[ "$6" != "uniref90" ]])
then
  mat="$6"
fi 

#Check the mode of operation: quiet/show
if [[ "$6" == "quiet" ]] || [[ "$7" == "quiet" ]] 
then
  mode="-q"
fi



#localizing fragments
locateFragment.pl -a  $1 -f $2 $mode -bdb $db
locateFragment.pl -a  $3 -f $4 $mode -bdb $db

#Aligning fragments and full sequences
alignSeqsFiles.pl -q $1_frag.faa -ql $1_frag -s $3_frag.faa -sl $3_frag -e 0.1 -c 20 -cc X -m $mat
alignSeqsFiles.pl -q $1.faa -ql $1 -s $3.faa -sl $3 -e 0.1 -c 5 -cc X -m $mat


#Open html reports if apropriate
if [[ $mode != "-q" ]]
then
  open ssearch*/*.html
fi
