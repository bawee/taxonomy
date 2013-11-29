#!/opt/local/bin/perl

# This Script is used together with CDS.IPR.blastp_parser.py
# Blast results that are used as imput for this program can be filtered (before running blast)  using a list of GIs ( retrieve_child_gis.pl )


use warnings;
use strict;

my $output;
my $Nodes_file = "nodes.dmp";
#my $Names_file = "names.dmp";
my $Names_file = "scientific_names.dmp"; #to only display scientific names, use this instead of the previous line.
my $GI_to_TAXID_file = "gi_taxid_prot.dmp";
my $TaxID;
my $refTaxID;

my @Taxonomy;
my @queryTaxIDHierachy;
my @refTaxonomy;
my @refTaxIDHierachy;
my %refTaxIDHash;

my $lastCommonGroup;

my $refGI = $ARGV[0]; #query GI number (i.e. GI number of BLAST hit)
my $queryGI = $ARGV[1]; #reference GI number (i.e. reference organism)




unless (scalar @ARGV == 2)	{
print <<USAGE and die;
	
USAGE: traverse_taxonomy.pl ReferenceGI QueryGI 

This script reports up the first common taxonomical group that is shared between a reference and query protein.
This script requires the files nodes.dmp, names.dmp and gi_taxid_prot.dmp in the working directory.
These files can be downloaded at ftp://ftp.ncbi.nih.gov/pub/taxonomy. nodes.dmp and names.dmp are found in taxdump.tar.gz 

USAGE
	
}



#my $refTaxID = get_name_or_taxid("$refGI", "$GI_to_TAXID_file", "2");

#THIS PART IS FOR THE REFERENCE

# Get the taxid corresponding to the GI number
$refTaxID = get_name_or_taxid("$refGI", "$GI_to_TAXID_file", "2");
chomp $refTaxID;

die "ERROR: Reference TaxID does not exist in the Protein databse. Could be nucleotide?\n" if $refTaxID eq ""; #checks for non-existant TaxIDs.

while ( $refTaxID > 1 )	{


#for (my $i = 1; $i < 10; $i++)	{
	
	#print $TaxID; #testing line
	# Obtain the scientific name corresponding to a taxid
	my $name = get_name_or_taxid("$refTaxID", "$Names_file", "3");
	chomp $name;
	#print $name;

	# Obtain the parent taxa Taxid
	my $parent = get_name_or_taxid("$refTaxID", "$Nodes_file", "3");
	chomp $parent;
	#print $parent;
	
	# Obtain the level name
	my $level = get_name_or_taxid("$refTaxID", "$Nodes_file", "5");
	chomp $level;
	#print $parent;
	
	# Build the taxonomy path
	unshift (@refTaxonomy, $name);
	unshift (@refTaxIDHierachy, "$refTaxID:$level:$name");
	$refTaxIDHash{ "$refTaxID" } = "$level";
	
	$refTaxID = $parent;
	
	
}

#print "For Reference: $refGI, the taxonomy is\n";
#print join("\n", @refTaxonomy);
#print join("\n", @refTaxIDHierachy), "\n";
#print "$Taxonomy[0]";


#print "\n\n-------------------------------------------------------------------------------------\n\n";


#THIS PART IS FOR THE QUERY
# Get the taxid corresponding to the GI number
$TaxID = get_name_or_taxid("$queryGI", "$GI_to_TAXID_file", "2");
chomp $TaxID;

die "ERROR: Query TaxID does not exist in the Protein databse. Could be nucleotide?\n" if $TaxID eq ""; #checks for non-existant TaxIDs.



while ( $TaxID > 1 )	{
#unless (exists($refTaxIDHash{$TaxID}))	{

#for (my $i = 1; $i < 10; $i++)	{

if (exists($refTaxIDHash{$TaxID}))	{ #checks if Query TaxID exists in the Taxonomy of the reference
	
	# Obtain the level name
	my $level = get_name_or_taxid("$TaxID", "$Nodes_file", "5");
	chomp $level;
	#print $parent;
	
	$lastCommonGroup = $level;
	
	last;
}

else	{
	
	#print $TaxID; #testing line
	# Obtain the scientific name corresponding to a taxid
	my $name = get_name_or_taxid("$TaxID", "$Names_file", "3");
	chomp $name;
	#print $name;

	# Obtain the parent taxa Taxid
	my $parent = get_name_or_taxid("$TaxID", "$Nodes_file", "3");
	chomp $parent;
	#print $parent;
	
	# Obtain the level name
	my $level = get_name_or_taxid("$TaxID", "$Nodes_file", "5");
	chomp $level;
	#print $parent;
	
	# Build the taxonomy path
	unshift (@Taxonomy, $name);
	unshift (@queryTaxIDHierachy, "$TaxID:$level:$name");
	$TaxID = $parent;
	
	
}
}



#print "For Query: $queryGI, the taxonomy is\n";
#print join("\n", @Taxonomy), "\n";
#print join("\n", @queryTaxIDHierachy), "\n";
#print "$Taxonomy[0]";
#print "\n\nThe query is related to the reference at the : $lastCommonGroup level \n";
print $lastCommonGroup; #This line needed by the colouring script

sub get_name_or_taxid # Obtain the name corresponding to a taxid or the taxid of the parent taxa
{
	my($a, $b, $c);
	($a, $b, $c) = @_;
	`grep --max-count=1 "^$a" "$b" | cut -f"$c"`; #this grep pulls out the corresponding taxid for the supplied GI
}



exit



