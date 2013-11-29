TaxIDermy
-----------------
a set of scripts compiled by Bryan Wee (1st Aug 2012)



Updated: 28th Nov 2012


Requirements
==============

BioPerl
Biopython
local copy of nr

1x .faa protein sequence multi fasta file from the genome that you want to analyse with locus tags
1x gbk file of the genome you want to analyse

1x gi_taxid_prot.dmp *this file is huge >1GB*
1x nodes.dmp
1x names.dmp
1x scientific_names.dmp *this file is derived from names.dmp. See Note 3*

Note:
	* The latest version of these files can be downloaded at ftp://ftp.ncbi.nih.gov/pub/taxonomy. nodes.dmp and names.dmp are found within taxdump.tar.gz
	* scientific_names.dmp is made using: `grep "scientific" names.dmp > scientific_names.dmp`
	* make sure the latest version of nodes.dmp, names.dmp and gi_taxid_prot.dmp are available so blast hits from the latest version of NR can be processed.

Instructions
===============

1. The first step is to create a list of gi numbers that correspond to sequences in NCBI that belong to the taxons that you do not want to include in the results.  This can be done by using the script retrieve_child_gis.pl "<Node to filter>"  and using the generated file <Taxon>.gi_list in the blastp step in step 2::

	$ retrieve_child_gis.pl "Pseudomonas"
	# or 
	$ retrieve_child_gis.pl "Pseudomonas aeruginosa PA14"

This script will output the gi numbers of associated with all taxies below the node "Pseudomonas" into a new file: <Taxon>.gi_list. The first line of the STD_OUT of this script will be <TaxonID> <Name of taxon> <Taxon level>. The rest of the STD_OUT of this script will be the list of Taxons that are below the Pseudomonas node.


2. Do a blastp using your the protein multi-fasta from the region you would like to colour against nr. Run blastp with the flag -negative_gilist <Pseudomonas>.gi_list::

	$ blastp  -query Pseudomonas_PA14.faa -db nr -outfmt 6 -out  Pseudomonas.v.nr.blastp.tab -num_threads 8 -negative_gilist Pseudomonas.gi_list


Note: nr must be in your path. if not, point the -db argument to nr on your computer


3. To create the final output, the following script needs to be run::

	$ python CDS.IPR.blastp_parser.py <blastp output> <Genome genbank> <TaxID of reference> <output.file.with.colour>
	$ python CDS.IPR.blastp_parser.py Pseudomonas.v.nr.blastp.tab Pseudomonas_PA14.gbk 446 > bryanscript.colour.out.tab

Note: this script calls traverse_taxonomy.pl which must be in the working directory
