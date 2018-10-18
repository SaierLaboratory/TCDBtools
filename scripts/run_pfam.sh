#!/bin/bash

#==========================================================================
# Arguments:
#   1. File with input sequences in fasta format
#   2. output file name

hmmscan --cpu 4 --noali --cut_ga -o /dev/null --domtblout $2 /ResearchData/pfam/pfamdb/Pfam-A.hmm $1
