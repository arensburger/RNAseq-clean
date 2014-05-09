#!/usr/bin/perl

# May 2014 Takes fastq library as input and returns quality statitics using fastx

use strict;
require File::Temp;
use File::Temp ();
use File::Basename;
use Getopt::Long;
use File::Path;

my $filename;
my $outputname; 

#####read and check the inputs
GetOptions(
	'i:s'   => \$filename,
	'o:s'	=> \$outputname
);

unless ($filename) {
	die ("usage: perl library_quality -i <REQUIRED: FASTQ file> -o <OPTIONAL output name>\n");
}
unless ($outputname) {
	$outputname = basename($filename, (".fq", ".fastq"));
}

# generate the stats
my $ts = File::Temp->new( UNLINK => 1, SUFFIX => '.txt' ); # temporary file
my $boxfilename = $outputname . "-quality_boxplot.png";
my $nucfilename = $outputname . "-nuc_dist.png";
`fastx_quality_stats -i $filename -o $ts -Q33`;
`fastq_quality_boxplot_graph.sh -i $ts -o "$boxfilename" -t $outputname`;
`fastx_nucleotide_distribution_graph.sh -i $ts -o "$nucfilename" -t $outputname`;
print "done";
exit;
