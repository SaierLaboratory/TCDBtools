#!/bin/bash

#Extract the clan info for any pfam accession
zgrep "$1" /ResearchData/pfam/download/Pfam-A.clans.tsv.gz

