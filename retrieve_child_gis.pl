#!/opt/local/bin/perl

# requires nodes.dmp and names.dmp in the working directory

use warnings;
use strict;

use Bio::DB::Taxonomy;

my $idx_dir = './tmp';

my %hash_of_child_taxids;
 
system "mkdir tmp";

my $species = $ARGV[0] or die "USAGE: xxx.pl <species>"; #user defines species name in quotes e.g. "Pseudomonas aeruginosa"


my ($nodefile,$namesfile) = ('nodes.dmp','names.dmp');
my $db = new Bio::DB::Taxonomy(-source    => 'flatfile',
                               -nodesfile => $nodefile,
                               -namesfile => $namesfile,
                               -directory => $idx_dir);

my $taxon = $db->get_taxon(-name => $species);

my $node = $db->get_Taxonomy_Node(-taxonid => $taxon->id); #this pulls out the taxid from the taxon specefied in the line above
print $node->id, " ", $node->scientific_name, " ", $node->rank, "\n";
get_descendents($node->id);


#print join("\n", @array_of_child_taxids);

my $index;

open FILE, "gi_taxid_prot.dmp" or die $!;
open WRITE, ">>$species.gi_list" or die $!;

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

exit;




sub get_descendents #This subroutine recursively searches for taxids below the given taxid until no taxids are found.
{


my $taxonid = $_[0];
my $node = $db->get_Taxonomy_Node(-taxonid => $taxonid);
my @extant_children = $db->each_Descendent($node);

if (@extant_children)	#checks if @extant_children is filled
	{  
	
		foreach my $child ( @extant_children )
		{
			get_descendents($child->id);
			#print "id is ", $child->id, "\n"; # NCBI taxa id
    		#print "rank is ", $child->rank, "\n"; # e.g. species
    		print $child->scientific_name, "\n"; # scientific name
			$hash_of_child_taxids{$child->id} = $child->scientific_name;
			#push(@array_of_child_taxids, $child->id); #previously used an array to store taxids. Updated to use hash for faster access
		}
	}
else	
	{
		return "";
	}

}
