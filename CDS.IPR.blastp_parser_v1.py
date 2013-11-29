#!/usr/bin/env python2

####################################
# Simon Belluzzo		
# s.belluzzo@student.unimelb.edu.au
#
# MODIFIED BY -
#   Mitchell Stanton-Cook 07/2011
#   To handle " and / chars in 
#   output (breaks artemis)
#   
#   Mitchell Stanton-Cook & Bryan Wee 24/08/2011
#   To parse the output of the web-based Batch CD-Search (format: full, with domain definition lines)
#   Major differences: Interpro "Confidence" value is replaced with CDSearch "incomplete" value
#
#   Bryan Wee 09/08/2012
#   CDS.IPR.blastp_parser.py To parse blastp output and display hits in genbank format.
#   
#   Bryan Wee 28/11/2012
#   1. Receives TaxID of reference instead of Gi number to feed into traverse taxonomy script
#   2. Receives the name of the taxon level (e.g. species Legionell_pneumophila) and prints it to %level
#
# TESTED: By Bryan Wee
# Usage: CDsearch_parser.py <CD search output file> <Corresponding genbank file> <CDS/IPR>
# 
#
####################################

import sys
import csv
import re
import optparse
from Bio import SeqIO
import subprocess

def read_iprscan(iprscanraw):
    hits = {} #initiates a dictionary
    hitsfile = csv.reader(iprscanraw, delimiter='\t')
    seen = []
    for hit in hitsfile:
        if hit[0] not in seen:
            seen.append(hit[0])
            sseqid = []	
            sseqid = hit[1].split('|')
            giNumber = sseqid[1]
            genbankAcc = sseqid[3]
            if hit[0] not in hits:
                hits[hit[0]] = [] #initiates array aka list in python. B.W:key seems to be name of database. 
            if len(hit) == 12: hits[hit[0]].append({'locus':hit[0],'acc':giNumber,'desc':genbankAcc,'start':hit[6],'end':hit[7],'evalue':hit[10],'score':hit[11]}) #Modified by BW to parse blastp tab output
            else: sys.stderr.write("Line does not have 12 columns.\n") #hits[hit[3]].append({'locus':hit[0],'acc':hit[4],'desc':hit[5],'start':hit[6],'end':hit[7],'score':hit[8],'conf':hit[9],'ipr':hit[11],'iprdesc':hit[12],'go':hit[13]})
    return hits
    
def read_genome_feat(genomefile):
    CDS = {}
    genome = SeqIO.read(genomefile,"genbank")
    for feature in genome.features:
        if feature.type == "CDS":
            #if not feature.strand: feature.strand = 1
            CDS[feature.qualifiers['locus_tag'][0]] = ({'strand' : feature.strand, 'start' : feature.location.nofuzzy_start, 'end' : feature.location.nofuzzy_end})
    return CDS


def write_out_tab(results,outfile):
    print "out"

def read_args():
    argparse = optparse.OptionParser()
    argparse.add_option("-g","--genome",dest="genbank",
                        help="GenBank or EMBL file to map hits to",metavar="FILE")
    argparse.add_option("-i","--input",dest="hits",
                        help="File of hits in tab format",metavar="FILE")
    argparse.parse_args()

def main():

    if len(sys.argv) != 4: #print USAGE message if no argv are given
        print """
This script parses the output of NCBI's blastp to a format that can be then read as an entry into Artemis.

Usage: CDS.IPR.blastp_parser.py <blastp_tab_output> <Corresponding genbank file>
 
 """ 
        sys.exit(1) 

    colours = {'no_rank':'201 31 22',
                'superkingdom':'217 88 14',
                'phylum':'252 207 3',
                'class':'206 213 0',
                'order':'105 176 17',
                'family':'8 131 67',
                'genus':'31 154 215',
                'species group':'17 50 121', 
                'species':'49 12 90',
                'species subgroup':'130 0 96',
                'subspecies':'197 0 89',
                '':'196 35 43'
            }

    hits = read_iprscan(open(sys.argv[1],'rU'))
    CDS = read_genome_feat(open(sys.argv[2],'rU'))
    ref_taxid = sys.argv[3]
    qlead = " "*21 # leading string to align feature qualifiers (adds spaces)
    for db in hits.keys(): # iterate through hits
        if db == "Seg": continue # skip Seg hits. Will make list from input to select dbs to use [BW: not applicable to CDSearch]
#        print db
        for hit in hits[db]:
            if hit['locus'] not in CDS: # check that gene the hit is against is present in genome file
                sys.stderr.write("%s not present in genome reference.\n" % hit[0]) 
                continue
            else: # actual working part
                if CDS[hit['locus']]['strand'] == -1: # if gene is on the -ve strand
                    hitend = int(CDS[hit['locus']]['end']) - int(hit['start'])*3   # swapped compared to normal,
                    hitstart = int(CDS[hit['locus']]['end']) - int(hit['end'])*3-2 # due to feature orientation
                    print "     misc_feature    complement(%i..%i)" % (hitstart, hitend)
                else: # all other cases (positive strand, not specified)
                    hitstart = int(CDS[hit['locus']]['start']) + int(hit['start'])*3-2
                    hitend = int(CDS[hit['locus']]['start']) + int(hit['end'])*3
                    print "     misc_feature    %i..%i" % (hitstart, hitend)
              
                print '%s/locus_tag="%s"' % (qlead,hit['locus'])
              
              #Here are the changes to handle / and " characters in iprdesc and desc - next 6 lines by MJSC & BW                    
                #this looks for and replaces illegal characters in iprdesc
                #fixed = hit['iprdesc'].replace('/', ' OR ')
                #fixed = fixed.replace('"', ' ')
                #hit['iprdesc'] = fixed
                #this looks for and replaces illegal characters in desc
                fixed = hit['desc'].replace('/', ' OR ')
                fixed = fixed.replace('"', ' ')
                hit['desc'] = fixed
                     

                print '%s/note="%s hit to %s, %s, evalue %s, bitscore %s"' % (qlead, db, hit['acc'], hit['desc'], hit['evalue'], hit['score'])
                
                traverse_call = subprocess.Popen(["perl traverse_taxonomy_v1.pl "+ref_taxid+" "+hit['acc']], shell=True, stdout=subprocess.PIPE)
                results_string = traverse_call.communicate()[0]
                results = results_string.split() #split the output of traverse_taxonomy.pl to 1.last commmon level and 2.taxa of level.
                #print results
                print '%s/colour=%s' % (qlead,(colours[results[0]])) #checks the colour from the colour dictionary above
                print '%s/common_node=%s' % (qlead, results[0]) #prints out the taxonomic level that the query has in common with the ref
                print '%s/origin="%s, from the same %s"' % (qlead, results[2], results[1]) #prints out the name of the taxonomic level that the query has in common with the ref




if __name__ == "__main__":
    main()
    
"""  
#Reference for the outputs of interproscan and blastp:

#qseqid   sseqid 					     pident length mismatch gapopen qstart qend sstart send   evalue bitscore
#CP1     gi|37955662|gb|AAP22501.1|      100.00  294     0       0       1    294     1   294     0.0      594
#  0              1                         2     3      4       5       6     7      8    9       10       11 


Locus            checksum       len   db         acc          desc     start   end     evalue              status  daterun  iprlookup
LPC_1272	006360BB50EF50B8	264	superfamily	SSF52096	SSF52096	1	   260	9.09992432541739E-74	T	20-Jul-2011	NULL	NULL
 0               1               2    3          4             5        6       7        8                  9        10     11       12     13


NF00181542	is the id of the input sequence.
27A9BBAC0587AB84	 is the crc64 (checksum) of the protein sequence (supposed to be unique).
272	 is the length of the sequence (in AA).
HMMPIR	is the anaysis method launched.
PIRSF001424	is the database members entry for this match.
Prephenate dehydratase	is the database member description for the entry.
1	is the start of the domain match.
270	is the end of the domain match.
6.5e-141	is the evalue of the match (reported by member database method).
T	is the status of the match (T: true, ?: unknown).
06-Aug-2005	is the date of the run.
IPR008237	is the corresponding InterPro entry (if iprlookup requested by the user).
Prephenate dehydratase with ACT region	is the description of the InterPro entry.
Molecular Function:prephenate dehydratase activity (GO:0004664)	is the GO (gene ontology) description for the InterPro entry.

"""
