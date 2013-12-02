#!/opt/local/bin/perl
#Written by Bryan Wee 1 August 2012
#updated by BW 2 Dec 2013

#This script takes a query phrase e.g. "Pseudomonas" or "Pseudomonas aeruginosa" or "Pseudomonas aeruginosa UCBP-PA14" and returns all GI numbers associated with ALL taxa beneath the query taxa.

#Requirements:
# requires an updated copy of nodes.dmp and names.dmp in the working directory

use warnings;
use strict;
use Bio::DB::Taxonomy;

#USAGE:
my $species = $ARGV[0] or die "USAGE: retrieve_child_gis.pl <species>\n\nWhere <species> is a query string such as species name or genus or family etc.\n\n"; #user defines species name in quotes e.g. "Pseudomonas aeruginosa"

my %hash_of_child_taxids; #initialize hash that stores all the taxids

my $idx_dir = './tmp'; #
system "mkdir tmp"; #create directory for the local taxonomy db

#this part initializes the database of the nodes and names used by Bio::DB::Taxonomy
my ($nodefile,$namesfile) = ('nodes.dmp','names.dmp');
my $db = new Bio::DB::Taxonomy(-source    => 'flatfile',
                               -nodesfile => $nodefile,
                               -namesfile => $namesfile,
                               -directory => $idx_dir);

my $taxon = $db->get_taxon(-name => $species) or die "\nInvalid query, please try again or use the scientific name...\n\n"; #retrieve taxon object from the local database using user query phrase. Fails if query is invalid

my $node = $db->get_Taxonomy_Node(-taxonid => $taxon->id); #this line pulls out the TaxID of the taxon object from the line above

$hash_of_child_taxids{$taxon->id} = $taxon->scientific_name; #this puts in the parent node/taxid into the hash. Previously this step was skipped resulting in GI numbers tied to the parent taxid to be missed. Added by BW 2.12.2013

print "The query phrase matches: TaxID ", $node->id, " (", $node->scientific_name, "), ", $node->rank, "\n\n"; #Prints out the taxID, scientific name and rank of the query

print "The taxa/taxon under the query taxon is/are:\n";
get_descendents($node->id); #This runs the main recursive subroutine define below

my $index;

$species =~ s/\s/_/g;
open FILE, "gi_taxid_prot.dmp" or die $!;
open WRITE, ">$species.gi_list" or die $!;
print "\nCollecting GI numbers.... (will take some time)\n";

while (my $line = <FILE>)	{ #for each line in gi_taxid_prot.dmp

	chomp $line;
	my @elements = split('\t', $line);
	my $taxid_gi2taxid = $elements[1];
	my $gi_gi2taxid = $elements[0];
	if (exists $hash_of_child_taxids{$taxid_gi2taxid}) #checks if exists in hash of child taxids
	{
		print WRITE $gi_gi2taxid, "\n";
		#`awk '\$2==($child){print \$1}' gi_taxid_prot.dmp >> $species.test.gi_list`; #used awk previously, took too long
	}
	else
	{
		next;
	}
}

close FILE;
close WRITE;

print "\nList of GI numbers associated with all taxa under $species have been printed to $species.gi_list\n";

exit;


sub get_descendents#This subroutine recursively searches for taxids below the given taxid until no further taxids are found.
{

my $taxonid = $_[0]; #receive the id of parent
my $node = $db->get_Taxonomy_Node(-taxonid => $taxonid); #get the node object from the local taxonomy db in ./tmp
my @extant_children = $db->each_Descendent($node); # also possible to use the method get_all_Descendents

if (@extant_children)	#checks if @extant_children is filled
	{  
	
		foreach my $child ( @extant_children )
		{
			get_descendents($child->id);
			#print "id is ", $child->id, "\n"; # NCBI taxa id
    		#print "rank is ", $child->rank, "\n"; # e.g. species
    		print $child->scientific_name, "\n"; # scientific name
			$hash_of_child_taxids{$child->id} = $child->scientific_name; #creates a hash with format taxid => scientific_name
			#push(@array_of_child_taxids, $child->id); #previously used an array to store taxids. Updated to use hash for faster access
		}
	}
else	
	{
		return "";
	}

}
