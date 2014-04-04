#!/opt/local/bin/perl

# This Script is used together with taxIDermy.py
# Blast results that are used as imput for this program can be filtered (before running blast)  using a list of GIs ( retrieve_child_gis.pl )


#update 2013 Nov 28 (Bryan Wee):
#	1. receives reference TaxID instead of GI as an argument
#	2. prints out name of last common taxon name as well as the level (e.g. "species Legionella_pneumophila_130b").
#	3. Tidied up script.


use warnings;
use strict;

my $output;
my $Nodes_file = "nodes.dmp";
#my $Names_file = "names.dmp";
my $Names_file = "scientific_names.dmp"; #to only display scientific names, use this instead of the previous line.
my $GI_to_TAXID_file = "gi_taxid_prot.dmp";
my $TaxID;
my $queryTaxon; 

my @Taxonomy;
my @queryTaxIDHierachy;
my @refTaxonomy;
my @refTaxIDHierachy;
my %refTaxIDHash;

my $lastCommonGroup;

my $refTaxID = $ARGV[0]; #query TaxID number (i.e. TaxID of organism that is used as reference.)
my $queryGI = $ARGV[1]; #query GI number (i.e. GI number from blast hit)


#open LOG, ">>logfile.txt" or die$!; #uncomment to generate logfile for testing purposes. Also uncomment additional lines that start with "print LOG" below


unless (scalar @ARGV == 2)	{
print <<USAGE and die;
	
USAGE: traverse_taxonomy.pl ReferenceGI QueryGI 

This script reports up the first common taxonomical group that is shared between a reference and query protein.
This script requires the files nodes.dmp, names.dmp and gi_taxid_prot.dmp in the working directory.
These files can be downloaded at ftp://ftp.ncbi.nih.gov/pub/taxonomy. nodes.dmp and names.dmp are found in taxdump.tar.gz 

USAGE
	
}

#PART1 : THIS PART pulls out the taxonomical hierarchy of the reference

# Get the taxid corresponding to the GI number
#$refTaxID = get_name_or_taxid("$refGI", "$GI_to_TAXID_file", "2"); #this line is for the GI of the reference instead of the taxid.
chomp $refTaxID;

if ($refTaxID eq "") {
	print "No_result No_result No_result"; #checks for non-existant TaxIDs.
	die "ERROR: Reference TaxID does not exist in the Protein databse. Could be nucleotide?\n";
}


while ( $refTaxID > 1 )	{ #keeps going until it reaches the root level i.e. TaxID = 1

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
	
	$refTaxID = $parent; #assigns the parent TaxID for the next iteration.

        if ($refTaxID == 1)        {
        unshift (@refTaxonomy, "root");
        unshift (@refTaxIDHierachy, "$refTaxID:root:root");	
	$refTaxIDHash{ "1" } = "root";
	}
	
}

#Extra parts to print out the entire taxonomical hierarchy stored in an array for testing or viewing

#print "For Reference: $refGI, the taxonomy is\n";
#print join("\n", @refTaxonomy);
#print LOG join("\n", @refTaxIDHierachy), "\n";
#print "$Taxonomy[0]";


#print "\n\n-------------------------------------------------------------------------------------\n\n";


#PART2 : HIS PART pulls out the taxonomical hierarchy of the query (blast hit)

# Get the taxid corresponding to the GI number of the blast hit specified in $ARGV[1]
$TaxID = get_name_or_taxid("$queryGI", "$GI_to_TAXID_file", "2");
chomp $TaxID;

if ($TaxID eq "") {
	print "no_rank No_result No_result"; #checks for non-existant TaxIDs.
	die "ERROR: Query TaxID does not exist in the Protein databse. Could be nucleotide?\n";
}

#die "ERROR: Query TaxID does not exist in the Protein databse. Could be nucleotide?\n" if $TaxID eq ""; #checks for non-existant TaxIDs.

$queryTaxon = get_name_or_taxid("$TaxID", "$Names_file", "3");

while ( $TaxID > 1 )	{ #keeps going until it reaches the root level i.e. TaxID = 1

#if (exists($refTaxIDHash{$TaxID}))	{ #checks if Query TaxID exists in the Taxonomy of the reference
	
	# Obtain the level name
#	my $level = get_name_or_taxid("$TaxID", "$Nodes_file", "5");
#	chomp $level;
	#print $parent;
#	$Taxonomy[0] = $level; #this line catches instances when the GI number is the same as $lastCommonGroup.
#	$lastCommonGroup = $level;
#	
#	last;
#}

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
	$TaxID = $parent; #assigns the parent TaxID for the next iteration.
	
	if ($TaxID == 1)	{
	unshift (@Taxonomy, "root");
        unshift (@queryTaxIDHierachy, "$TaxID:root:root");
	}

	if (exists($refTaxIDHash{$TaxID}))      { #checks if Query TaxID exists in the Taxonomy of the reference

       		# Obtain the level name
       		my $level = get_name_or_taxid("$TaxID", "$Nodes_file", "5");
      	 	chomp $level;
       		#print $parent;
		$Taxonomy[0] = $level; #this line catches instances when the GI number is the same as $lastCommonGroup.
	       	$lastCommonGroup = $level;

        	last;
	}


}

#Extra parts to print out the entire taxonomical hierarchy stored in an array for testing or viewing

#print LOG "For Query: $queryGI, the taxonomy is\n";
#print LOG join("\n", @Taxonomy), "\n";
#print LOG join("\n", @queryTaxIDHierachy), "\n";
#print "$Taxonomy[0]";
#print LOG "\n\nThe query is related to the reference at the : $lastCommonGroup level \n";

$queryTaxon =~ s/\s/_/g;
$Taxonomy[0] =~ s/\s/_/g; #replaces spaces in "no rank" with an underscore
$lastCommonGroup =~ s/\s/_/g; #replaces spaces in taxon tame with underscores
print $lastCommonGroup , " ", $Taxonomy[0], " ", $queryTaxon; #This prints the STDOUT that is needed by taxIDermy.py

sub get_name_or_taxid # Obtain the name corresponding to a taxid or the taxid of the parent taxa
{
	my($a, $b, $c);
	($a, $b, $c) = @_;
	`grep --max-count=1 "^$a" "$b" | cut -f"$c"`; #this grep pulls out the corresponding taxid for the supplied GI
}


#close LOG; #closes logfile used for testing purposes
exit



