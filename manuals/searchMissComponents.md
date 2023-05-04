# Documentation for script: _searchMissComponents.py_


## Summary
When characterizing a multicomponent system in a (meta)genome, often some components in TCDB will not have good 
matches with proteins in the (meta)genome under study. This script takes functional keywords associated with the missing 
component and downloads from NCBI all non-redundant proteins containing those keywords in their annotations. Sequences 
are then BLASTed against the query (meta)genome and for all significant hits hydropathy plots and Pfam domains are
shown to help the user determine whether a given missing component was found.

## Contributors  
Yichi Zhang and Arturo Medrano-Soto


## Dependencies
The following Python module needs to be available to the script: 

1. **Python 2.X**  
You can download the Python 2 from the [official website](https://www.python.org/). The following
modules are required:  
  a) Biopython 1.75  
  b) pandas 0.24.2  
  c) requests 2.24.0  
  d) formatQuodDomain (included in the BioVx distribution)  
  e) hmmscanParser (included in the BioVx distribution)  

2. **_blast+ 2.6.0 to 2.10.0_**  
Other versions of blast may require minor adaptations. Visit the
[download site](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download).  

3. **_HHMER 3.2.1_**  
You can download this suite of programs from its [official site](http://hmmer.org/).  

4. **_quod.py_**  
This program in part of this [distribution](https://gitlab.com/khendarg/hvordan/blob/master/docs/quod.md)..  

5. **Proteome sequence file**   
File with all sequences in FASTA format of the (meta)genome under study.  


## Command line options
The following options are available:

    usage: searchMissComponents.py [-h] -k <keywords string> [-l] -g <string
                               input_file> [-e <evalue>]
                               [-cf <a number between 0-100>] [-w <x,q,s,b>]
                               [-ef <evalue threashold for the result>]
                               [-dc <a number between 0-100>]
                               [-al <minimum alignment kept by the blast result>]
                               [-r <sequence identity threashold>]
                               [-d DISPLAY] [-out <string output_directory>]
                               [-o]
    Argument description:
    -h, --help            show this help message and exit
    
    -k <keywords string>, --keyword <keywords string>
                        MANDATORY. Descriptions of the missing components for
                        searching on the NCBI website. Logical expressions
                        like "AND","OR" can be used for detailed description.
                        
    -l, --local         allow user to choose if the sequences are retrieved
                        from a local NR blast database. The default setting is
                        to download from NCBI online. If the number of sequences
                        is much smaller than what were found online, it will
                        retrieve from the local database. The final result
                        will be the one with more sequences.
                        
    -g <string input_file>, --genome <string input_file>
                        MANDATORY. Name or location of a complete genome file
                        in FASTA format or FASTA compressed format(bz,gz).
                        
    -e <evalue>, --evalue <evalue>
                        the e-value threashold for blast. Default value is
                        1e-6. Users can lower the e-value if they found
                        nothing from blast.
                        
    -cf <a number between 0-100>, --coverage_filter <a number between 0-100>
                        percent coverage of the alignment of the smaller
                        protein, used for filtering blast results. default is
                        50
                        
    -w <x,q,s,b>, --which <x,q,s,b>
                        decide whether coverage filter will be applied on
                        query(q), subject(s), both(b) or either(x) one.
                        default is "x"
                        
    -ef <evalue threashold for the result>, --evalue_filter <evalue threashold for the result>
                        an additional filter for blast results. This filter
                        ensures that e-values of all blast results must be
                        equal to or smaller than the threshold. Default
                        threshold is 1e-6
                        
    -dc <a number between 0-100>, --domain_coverage <a number between 0-100>
                        a domain filter for pfam results. Default is 35
                        percent. Domains which coverages are over the
                        threshold will be considered as true hits.
                        
    -al <minimum alignment kept by the blast result>, --alignment_length <minimum alignment kept by the blast result>
                        a filter for blast results that keeps alignment length
                        (without gaps) greater than the threshold. Default
                        value is 30
                        
    -r <sequence identity threashold>
                        a global sequence identity threashold for subject
                        sequences. Default is 0.9. It removes redundant
                        sequences whose similarities are higher than 0.9
                        
    -d DISPLAY, --display DISPLAY
                        the minimum number of results displayed in the final
                        html output. Default is 100. if the number of blast
                        result is less than the threshold, it will display all
                        results
                        
    -out <string output_directory>
                        the name, or path and name of the output directory
                        containing all files generated by this program and the
                        final result in html format. Default name is
                        "MissingComponent_result" generated in the current
                        directory
                        
    -o, --overwrite     overwrite all output files. Disabled by default
