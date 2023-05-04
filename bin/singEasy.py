#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import os.path
import os
from Bio import SeqIO
from hmmscanParser import hmmParser
import argparse
import textwrap
import sys
header = ["tname", "tacc", "tlen", "qname","qacc", "qlen", "E-value", "score", "bias", "#", "of","c-Evalue", "i-Evalue", "score", "bias", "hfrom", "hto","afrom", "ato", "efrom", "eto", "acc", "description of target"]
dirpath = os.getcwd()


# In[2]:


def get_systems(tcdb):
    '''
    input a tcdb fasta file and return two lists. One for all single-component systems, one for multi-component systems
    '''
    mc = []
    sc = []
    systems = dict()
    with open(tcdb, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            tcid = record.id.split('-')[0]
            if tcid not in systems:
                systems[tcid] = 1
            else:
                systems[tcid] = systems[tcid] + 1

    for key in systems:
        if systems[key] == 1:
            sc.append(key)
        else:
            mc.append(key)

    return sc,mc

def common_domain(qacc, sacc, xgbhit_dir, dshare):
    comparison = 'ssearch_{}_vs_{}'.format(qacc,sacc)
    pfam_path = os.path.join(xgbhit_dir,qacc,comparison,'files','hmmscan.out')
    hmm = hmmParser(pfam_path)
    df=pd.DataFrame(hmm.matrix,columns=header)
    df=df[df.columns[:4]]
    domain = dict()
    for index, row in df.iterrows():
        if row['tacc'] not in domain:
            domain[row['tacc']] = [row['qname']]
        else:
            domain[row['tacc']].append(row['qname'])
    # at least one common domain
    dnum = 0
    for key in domain:
        if len(domain[key]) > 1:
            dnum = dnum + 1
    if dnum >= dshare:
        return True
    else:
        return False

def format_html(df,html_name,xgbhit_path):
    outfile = open(os.path.join(dirpath,html_name),'w')
    outfile.write('<html><table width="100%%" border="1"> <tr> <td>Query ID</td> <td>Hit ID</td><td>Hit TCID</td> <td>Hit Description</td> <td>Match Len</td> <td>e-Val</td> <td>% Identity</td> <td>Query Length</td> <td>Hit Length</td> <td>Query Coverage</td> <td>Hit Coverage</td> <td>Query TMS#</td><td>Hit TMS#</td> <td>TM-Overlap Score</td> <td>Family Abrv.</td><td>Predicted Substrate</td> <td> #</td> </tr><br>')
    for index,row in df.iterrows():
        outfile.write('<tr>')
        plot_path = os.path.join(xgbhit_path,row['#Query_id'],'ssearch_{}_vs_{}'.format(row['#Query_id'],row['Hit_xid']),'report.html')
        outfile.write('<td><a href={} target="_blank">{}</a></td>'.format(plot_path,row['#Query_id']))
        outfile.write('<td>{}</td>'.format(row['Hit_xid']))
        outfile.write('<td><a href="http://tcdb.org/search/result.php?tc={}">{}</a></td>'.format(row['Hit_tcid'],row['Hit_tcid']))
        outfile.write('<td>{}</td>'.format(row['Hit_desc']))
        outfile.write('<td>{}</td>'.format(row['Match_length']))
        outfile.write('<td>{}</td>'.format('{:.2e}'.format(row['e-value'])))
        outfile.write('<td>{}</td>'.format(row['%_identity']))
        outfile.write('<td>{}</td>'.format(row['Query_Length']))
        outfile.write('<td>{}</td>'.format(row['Hit_Length']))
        outfile.write('<td>{}</td>'.format(row['Query_Coverage']))
        outfile.write('<td>{}</td>'.format(row['Hit_Coverage']))
        outfile.write('<td>{}</td>'.format(row['Query_n_TMS']))
        outfile.write('<td>{}</td>'.format(row['Hit_n_TMS']))
        outfile.write('<td>{}</td>'.format(row['TM_Overlap_Score']))
        outfile.write('<td>{}</td>'.format(row['Family_Abrv']))
        outfile.write('<td>{}</td>'.format(row['Predicted_Substrate']))
        outfile.write('<td>{}</td>'.format(row['row_number']))
        outfile.write('</tr>')   


# In[3]:


def process_data(rawdata, tcdb_file, xgbhit_dir, evalue, hcov, qcov, dshare):
    df = pd.read_csv(rawdata,sep='\t')
    sc,mc = get_systems(tcdb_file)
    df_single = df[df['Hit_tcid'].isin(sc)]
    df_multi = df[df['Hit_tcid'].isin(mc)]
    # apply filter for all single-components systems
    # based on e-value
    df_single = df_single.astype({'e-value': 'float','Query_Coverage':'float','Hit_Coverage':'float'})
    df_single = df_single[(df_single['e-value'] <= float(evalue)) & (df_single['Query_Coverage'] >= float(qcov)) & (df_single['Hit_Coverage'] >= float(hcov))]
    # check domains
    keep_list = []
    keep_tcid = []
    for index,row in df_single.iterrows():
        #print(row['Hit_xid'])
        if common_domain(row['#Query_id'],row['Hit_xid'],xgbhit_dir,dshare):
            keep_list.append(index)
            keep_tcid.append(row['Hit_tcid'])
        if row['Hit_tcid'] in keep_tcid:
            keep_list.append(index)
    # get clean separate tables
    df_single_con = df_single[df_single['Hit_tcid'].isin(keep_tcid)]
    
    df_single_sus = df[df['Hit_tcid'].isin(sc)]
    df_single_sus.drop(list(df_single.index),inplace=True)
    # sort based on same tcid and e-value
    df_single_con = df_single_con.sort_values(by=['Hit_tcid','e-value'],ascending=True)
    df_single_sus = df_single_sus.sort_values(by=['Hit_tcid','e-value'],ascending=True)
    df_multi = df_multi.sort_values(by=['Hit_tcid','e-value'],ascending=True)
    # output a tsv sheet
    df_single.to_csv(os.path.join(dirpath,'confident_single.tsv'), sep='\t')
    df_single_sus.to_csv(os.path.join(dirpath,'suspects_single.tsv'), sep='\t')
    df_multi.to_csv(os.path.join(dirpath,'suspects_multi.tsv'), sep='\t')
    format_html(df_single_sus,'suspects_single.html',xgbhit_dir)
    format_html(df_single,'confident_single.html',xgbhit_dir)
    #return df_single, df_single_sus


# process_data('./results.tsv','./tcdb.faa','/ResearchData/Users/amedrano/RalfRabus/GBLAST/GCA_001854245.1_ASM185424v1/xgbhit','1e-10',70,70,1)

# In[ ]:


if __name__ == '__main__':
    parser = argparse.ArgumentParser(epilog=textwrap.dedent('''Description: 
    This program extracts all single-component transport systems from a raw GBlast output table. 
    Results are reported in two tab-delimited tables: 1) the most reliable assignments, and 2) less reliable assignments that need human validation.'''))
    # create a parser to handle options and arguments. Show help info if no args
    parser.add_argument( '-d', '--tabledata', type = str, dest = 'data', required = True, metavar = '<GBlast output file>', help = 'MANDATORY. Path to the output table in tsv format generated by the program GSAT. This table contains the top matches between the query genome and systems in TCDB.' )
    parser.add_argument( '-db', '--database', type = str, dest = 'db', required = True, metavar = '<sequence file>', help = 'MANDATORY. Path to the file in fasta format with All the protein content in TCDB as generated by the program extractFamily.pl. This file should reflect the TCDB version that will be used to analyze both single and multicomponet systems in the query genome.')
    parser.add_argument( '-ad', '--address', type = str, dest = 'ad', required = True, metavar = '<plots directory>', help = 'MANDATORY. Path to the root directory containing all hydropathy plots generated by GBlast.' )
    parser.add_argument( '-e', '--evalue', type = float, dest = 'evalue', default = '1e-10', help = 'The E-value cutoff for infering GBlast hits as homologs. By default E-value <= 1e-10.')
    parser.add_argument( '-scov', '--subjectcoverage', type = float, dest = 'hcov', default = 70.0, help = 'The alignment coverage cutoff for subject proteins. Default is 70 percent.')
    parser.add_argument( '-qcov', '--querycoverage', type = float, dest = 'qcov', default = 70.0, help = 'The alignment coverage cutoff for query proteins. Default is 70 percent.')
    parser.add_argument( '-ds', '--domainshared', type = int, dest = 'dshare', default = 1, help = 'The minimum number of shared Pfam domains between query and subject proteins. Default is 1.')
    
    args = parser.parse_args()
    if len(sys.argv) < 3:
        parser.print_help()
        sys.exit(1)
    
    process_data(args.data, args.db, args.ad, args.evalue, args.hcov, args.qcov, args.dshare)

