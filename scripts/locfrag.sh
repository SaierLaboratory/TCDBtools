
#Help
if [ "$1" = "-h" ]
then
  echo "Locate aligned fragments within 2 full proteins"
  echo "Arguments:"
  echo " 1. NCBI Accession of first protein 1"
  echo " 2. Aligned fragment of protein 1" 
  echo " 3. NCBI Accession of first protein 2"
  echo " 4. Aligned fragment of protein 2"
  echo " 5. (Optional) substitution matrix to use. (Defaul: BL50)"
  echo " 6. (Optional) Indicate whether plots will be shown (Values: show/quiet; Default: show)" 
  exit 1
fi



#Define the substitution matrix to work with
mat="BL50"
mode=""


#Identify the type of substitution matrix, if given
if [[ ! -z "$5" ]] && ([[ "$5" != "quiet" ]] && [[ "$5" != "show" ]])
then
  mat=$5
fi 


#Check the mode of operation: quiet/show
if [[ "$5" == "quiet" ]] || [[ "$6" == "quiet" ]] 
then
  mode="-q"
fi

locateFragment.pl -a  $1 -f $2 $mode 
locateFragment.pl -a  $3 -f $4 $mode
alignSeqsFiles.pl -q $1_frag.faa -ql $1_frag -s $3_frag.faa -sl $3_frag -e 0.1 -c 20 -cc X -m $mat
alignSeqsFiles.pl -q $1.faa -ql $1 -s $3.faa -sl $3 -e 0.1 -c 5 -cc X -m $mat


#Open html reports if apropriate
if [[ $mode != "-q" ]]
then
  open ssearch*/*.html
fi
