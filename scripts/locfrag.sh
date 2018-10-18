locateFragment.pl -a  $1 -f $2 
locateFragment.pl -a  $3 -f $4
alignSeqsFiles.pl -q $1_frag.faa -ql $1_frag -s $3_frag.faa -sl $3_frag -e 10 -c 20 -cc X
alignSeqsFiles.pl -q $1.faa -ql $1 -s $3.faa -sl $3 -e 10 -c 5 -cc X
open ssearch*/*.html
