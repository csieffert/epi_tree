#!/usr/bin/env perl
use strict;
use warnings;

use FindBin;
use lib $FindBin::Bin.'/../lib';
use Bio::Phylo::IO;
use Bio::Phylo::Forest::Tree;
use Bio::TreeIO;
use Text::CSV;
use Pod::Usage;
use Getopt::Long;
use File::Temp 'tempdir';
use Pod::Usage;
use Logger;

#=============================================================================
#Purpose:
#    Re-roots the phylogenetic tree with the indicated strain.
#
#Input: 
#    $input_taxa_tree -> Bio::Phylo::Forest::Tree phylogenetic tree to re-root
#    $newRootStrain -> String name of the strain to root with
#    $logger->reference to the Logger object.
#Output: 
#    A modified phylogenetic tree, re-rooted on the $newRootStrain
#=============================================================================
sub reRootTree
{
    my ($input_taxa_tree, $newRootStrain, $logger) = @_;

    my $newRoot = "'".$newRootStrain."'";
    foreach my $node ($input_taxa_tree->get_nodes()){ 
        if($node->id eq $newRoot)
        {
       	$input_taxa_tree->reroot_at_midpoint($node, 'ROOT NODE');
		$logger->log("The phylogenetic tree has been successfully re-rooted on strain: ".$node->id()."\n", 0);
		return;	
       }
    }   
    $logger->log("The requested strain".$newRoot."could not be found in the phylogenetic tree.\n", 0);
}

#==============================================================================
#Purpose:
#	Submodule to rearrange the entries in matrix.csv to match the new phylogenetic 
#ordering.
#Input: 
#	$input_taxa_tree-> Bio::Phylo::Forest::Tree phylogenetic tree to re-root
#	$inputMatrixFile-> The location of the input matrix.csv file.
#	$logger->reference to the Logger object.
#output: 
#   New matrix.csv file that matches the ordering of the re-rooted phylo tree
#==============================================================================
sub updateMatrixCsv
{
	my ($input_taxa_tree, $inputMatrixFile, $logger) = @_;
	
	#open a new file handle to print the output to
	open(my $revisedMatrixCsv, '>revisedMatrix.csv') or die "ERROR: Could not open the output file: $!";
	
	#open file handle for input matrix.csv file
	open(my $data, '<', $inputMatrixFile) or die "ERROR: Could not open '$inputMatrixFile' $!\n";
	
	my $csv = Text::CSV->new({ sep_char => '\t' });
	#hash the two-dimensional matrix as 'keyColumn:keyRow' : 'value' pairs to facilitate rearrangement 
	my %matrixHash = ();
	my @strainColumn=[];
	#parse the input matrix and retain the first row for strain names in @strainNames
	while (my $input = <$data>) {
		my @line = split(/\t/, $input);
		#add all of the strain names for each row
		if($line[0] eq 'strain'){
			@strainColumn = split(/ /, "@line");
		}
		else{
			my $strainRow = $line[0];
			my $index = 0;
			#hash the values in the matrix as 'keyColumn:keyRow' : 'value'
			foreach(@line){
				$matrixHash{"'".$strainColumn[$index]."'".':'."'".$strainRow."'"} = $_; 
				$index++; 
			}	
		}
	}
			
	#using reference tree, print a re-ordered matrix.csv file to the indicated file handle
	print $revisedMatrixCsv "strain\t";
	foreach($input_taxa_tree->get_nodes){
		if($_->is_Leaf()){
			print $revisedMatrixCsv $_->id."\t";
		}
	}
	print $revisedMatrixCsv "\n";
	foreach( $input_taxa_tree->get_nodes) {
		if($_->is_Leaf()){
    		my $rowNode = $_->id;
    		print $revisedMatrixCsv $_->id."\t";
    		foreach( $input_taxa_tree->get_nodes) {
    			if($_->is_Leaf()){
    				my $hashQuery = $rowNode.':'.$_->id;
    				my $hashResult = $matrixHash{$hashQuery};
    				print $revisedMatrixCsv $hashResult."\t";
    			}
    		}
    		print $revisedMatrixCsv "\n";	      
		}
    }
	close($revisedMatrixCsv);
	close($data); 
}

#==============================================================================
#Purpose:
#	Converts the branch lengths to an estimate of the total number of SNP 
#differences and renames the nodes to match the format: '[STRAIN][[BRANCH_LENGTH], [SNP_ESTIMATE]]' 
#Input:
#	$input_taxa_tree -> Bio::Phylo::Forest::Tree phylogenetic tree to re-root
#	$inputPhyFile -> The location of the input pseudoalign.phy file.
#	$logger->reference to the Logger object.
#==============================================================================
sub branchLengthToSNP
{
	my ($input_taxa_tree, $inputPhyFile, $logger) = @_;
		
	#parse the total SNP's in the tree from the input .phy file:
	open(my $inputPhy, "<", $inputPhyFile) or die "ERROR: Could not open input phylogeny.phy file.";
	my $input = <$inputPhy>;
	my @line = split(/\s/, $input);	
	
	my $treeTotalSNP = $line[2];
	print $treeTotalSNP."\n";
	my $internalNumber = 1;
	foreach my $node ( $input_taxa_tree->get_nodes ){
	  my $nodeBranchLength = $node->branch_length();
	  $nodeBranchLength = 0 if !defined $nodeBranchLength;
	  my $lengthToSNP = $nodeBranchLength*$treeTotalSNP;
	  $lengthToSNP = sprintf "%.2f", $lengthToSNP;
	  my $nodeName = $node->id();
      
      if($lengthToSNP<0.5){
      	$node->branch_length(0);
      }
      else{
      	$node->branch_length($lengthToSNP);
      }
      print $node->id() if $node->is_Leaf;
      $internalNumber++ if !$node->is_Leaf;
    }
    $logger->log("Internal node branches have been re-labelled to show total estimated SNP differences.\n", 0);
    close($inputPhy);
}

#==============================================================================
#Purpose:
#	Exponentially scales branch lengths on the input tree to allow the user to 
#	make the final output more easily human readable.
#Input:
#	$input_taxa_tree -> Bio::Phylo::Forest::Tree phylogenetic tree to re-root
#   $exponent -> the exponent to factor all branch lengths by
#	$logger -> Reference to a Logger object.
#Output:
#	matrix.csv file with resized branch lengths.
#==============================================================================
sub resizeTree
{
	my ($input_taxa_tree, $exponent, $logger) = @_;
	
	my $minimumBranchLength = 0.009;
	#zero out the branch values that are lower than threshold to maintain tree structure for near 0 branch lengths:
	foreach( $input_taxa_tree->get_nodes){
		if($_->get_branch_length() < $minimumBranchLength){
			$_->set_branch_length(0);
		}
	}
	
	$input_taxa_tree->exponentiate($exponent);
	$logger->log("Branch lengths have been resized with an exponential factor of: ".$exponent."\n", 0);
}

#==============================================================================
#Purpose:
#	Script to re-root a phylogenetic tree, order the tree in increasing/
#	decreasing order, convert branch lengths to total SNP's, and output a 
#	revised matrix.csv to match the phylogenetic tree.
#==============================================================================

my ($tmp_dir, $keep_tmp, $root_strain, $tree_order, $input_dir, $output_dir, $matrix_input, $input_phy, $help );

GetOptions(
	't|tmp-dir=s' => \$tmp_dir,
	'k|keep_tmp=s' => \$keep_tmp,
	'r|root=s' => \$root_strain,
	's|sort=s' => \$tree_order,
	'i|inputdir=s' => \$input_dir,
	'o|outdir=s' => \$output_dir,
	'm|matrix=s' => \$matrix_input,
	'p|phy=s' =>	\$input_phy,
	'h|help' => \$help
);
#handle documentation and help
if ($help){
	pod2usage(	-verbose => 99,
				-sections => "SYNOPSIS|OPTIONS|DESCRIPTION");
}

#check to ensure all required command line variables are present and set default values if applicable:
pod2usage(1) unless $tmp_dir && $input_dir && $output_dir && $matrix_input && $input_phy;
die "Error: No temp directory defined." if (not defined $tmp_dir);
$keep_tmp = 0 if (not defined $keep_tmp);

my $job_out = tempdir('rearrange_snp_matrixTempXXXX', CLEANUP => (not $keep_tmp), DIR => $tmp_dir) or die "Could not create temp directory";
my $logger = Logger->new($job_out);
	
#parse the newick format phylogeny generated by phyml into a Bio::Phylo::Forest::Tree object
my $taxa_file = $input_dir.'/pseudoalign.phy_phyml_tree.txt';
my $taxaIO = Bio::TreeIO->new(
	'-format'   => 'newick',
   	'-file'   => $taxa_file
);
my $taxa = $taxaIO->next_tree;
 	
#reroot the tree if requested by user:
reRootTree($taxa, $root_strain, $logger) if defined $root_strain;
 	
#parse the Bio::Tree::Tree object into a Bio::Phylo::Forest::Tree object to facilitate reordering
my $tree = Bio::Phylo::IO->parse(
   	'-string' => $taxa->as_text('newick'),
   	'-format' => 'newick'
)->first;
 	
#determine whether the tree should be re-sorted in decreasing or increasing order
if(defined $tree_order){
	if($tree_order eq "decreasing"){
		$tree->ladderize();
		$logger->log("The tree has been successful sorted in decreasing order.");
	}
	elsif($tree_order eq "increasing"){
		$tree->ladderize(1);
		$logger->log("The tree has been successful sorted in increasing order.");
	}
}
#create a matrix.csv file to reflect the changes made to the phylogenetic tree
#updateMatrixCsv($tree, $matrix_input, $logger) if (defined $root_strain || defined $tree_order);
#add the branch length to SNP labels to the tree nodes that show as: STRAIN[BRANCH_LENGTH, TOTAL_SNP's]
branchLengthToSNP($tree, $input_phy, $logger);
 	
#print the final newick formatted tree to a file that can be opened by any tree viewing program
open(my $treeout, '>phylogeneticTree.txt') or die "Could not open output file: $!";
print $treeout $tree->to_newick( -nodelabels => 1, -header => 1, -links => 1 );
close($treeout);

1;

=head1 NAME

rearrange_snp_matrix.pl - Script to re-root a phylogenetic tree, order the tree in increasing/decreasing order, convert branch lengths to total SNP's, and output a revised matrix.csv to match the phylogenetic tree. 

=head1 VERSION

This documentation refers to rearrange_snp_matrix.pl version 0.0.1.

=head1 SYNOPSIS

 rearrange_snp_matrix.pl -t tempDir -k keepTemp -r rootStrain -s sortOrder -i inputDataDirectory -o outputDir -m matrix.csv_input -p pseudoalign.phy_input

=head1 OPTIONS

=over

=item B<-t>, B<--tmp-dir> [required]

The directory location of the temporary directory where log messages are sent.

=item B<-k>, B<--keep-tmp> [optional]

Boolean flag indicating whether to keep or delete the temporary directory upon exiting the script.

=item B<-r>, B<--root> [optional] 

The output file.

=item B<-s>, B<--sort> [optional]

Flag, 'increasing' or 'decreasing', indicating the order in which to sort nodes in the output phylogenetic tree.

=item B<-i>, B<--inputDir> [required]

Directory containing the output files from phyml.

=item B<-o>, B<--outDir> [required]

The directory where all output files are located.

=item B<-p>, B<--phy> [required]

The location of the pseudoalign.phy file. 

=item B<-m>, B<--matrix> [required]

The location of the matrix.csv file.

=item B<-h>, B<--help>

To display help message

=back

=head1 DESCRIPTION

rearrange_snp_matrix steps:

1.  Takes a newick formatted phylogenetic tree file and reroots the tree on the strain indicated by --root.  The tree can then be sorted in increasing or decreasing order and branch lengths are converted and displayed as the total number of SNP's for each branch in the tree.  The new phylogenetic tree is output in newick format as a text file named phylogeneticTree.txt.

2.  A new matrix.csv file is generated that matches the ordering of the phylogenetic tree done in step 1 above.

=cut