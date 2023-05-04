# Documentation for script: _compareGblast.py_

## Summary
When using _GBlast_ to identify transport proteins in multiple genomes, a comparison
between the outputs of _GBlast_ may yield some insight into potential relationships between organisms.
This program compares multiple outputs of _GBlast_ and generates a tab-seperated table
summarizing the best hits from each _GBlast_ output for each TCDB protein.


## Contributors
Vasu Iddamsetty and Arturo Medrano-Soto.

## Dependencies
The program requires the following in the path to run properly:

1. **_Python 2.7+_**
The program was written to run with Python 2.7.14 but may work with more
recent versions of Python.


## Command Line Options
The following options are required for the program to run properly:

    positional arguments:
      files                 A space seperated list of paths to the relevant GBLAST
                            files.

    optional arguments:
       -h, --help            show this help message and exit
       -out OUTPUT, --output OUTPUT
                            REQUIRED. The path to the outfile.
       -id ID [ID ...]       REQUIRED. The space separated list of names that will
                            be used in the header to identify each genome.Should
                            be in the same order as the filenames.  
