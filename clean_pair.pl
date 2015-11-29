#!/usr/bin/perl
# Wed 21 Mar 2012 02:30:34 PM PDT automates initial library prep
# Fri 23 Mar 2012 10:17:03 AM PDT modified to work with a pair of libraries
# Wed 04 Jul 2012 08:07:18 AM PDT modified to work on hbar, removed Temp:seekable dependance and prevented making graphs because system lacks "gnuplot" library
# Fri 06 Jul 2012 10:57:13 AM PDT added a test to see if data is ok for chastity filter
# Sept/Oct 2012 modified the ribosome matching to accomodate spaces
# Fri 05 Oct 2012 12:59:25 PM PDT changed output to comma separated
# April 2014 now using SortMeRNA as ribosome removal program, updated Trimmomatic

use strict;
require File::Temp;
use File::Temp ();
#use File::Temp qw/ :seekable /;
use File::Basename;
use Getopt::Long;
use File::Path;

# CONSTANTS
# adapter trimming parameters
my $TRIMMOMATIC_PATH = "./Trimmomatic-0.32"; # java program location
my $SORTMELOC = "/home/peter/RNAseq-clean/sortmerna"; # must keep full path name for sortme program
my $MINLEN = 35;

# other parameters
my $readspair1; #filename of first pair fastq
my $readspair2; #filename of second pair fastq file
my $outputname; #base name for output files
my $outputdir; #outputdirectory
my $threads = `grep -c processor /proc/cpuinfo`; #number of threads to use
$threads =~ s/\s//g
;
#return date and time
sub datetime {
	use POSIX qw/strftime/;
	return (strftime('%D %T',localtime));
}

#####read and check the inputs
GetOptions(
	'1:s'   => \$readspair1,
	'2:s'	=> \$readspair2,
	'o:s'	=> \$outputname,
	'd:s'	=> \$outputdir,
	't:s'   => \$threads,
);

unless ($readspair1 and $readspair2) {
	die ("usage: perl clean_pair.pl -1 <FASTQ file pair 1> -2 <FASTQ file pair 2> -o <OPTIONAL output name> -d <OPIONAL output directory (default current directory)> -t <OPTIONAL: number of threads to use, default $threads>\n");
}
unless ($outputname) {
	$outputname = basename($readspair1, ".fastq"); #take the name from first member of the pair
}

#test if output directory has been specified, if not set output to current directoy
unless ($outputdir) {
	$outputdir = `pwd`; #set output to current directoy
	chomp $outputdir;
}

#test if output directory exists if not, create it
unless ( -e $outputdir )
{
	my $dir_test = mkdir($outputdir, 0777);
	unless ($dir_test) { 
		die "cannot create directory $outputdir\n";
	}
}

#create temporary files
my $paired_output = File::Temp->new( UNLINK => 1, SUFFIX => '.fastq' ); # temporary file with results pair1
my $unpaired_output = File::Temp->new( UNLINK => 1, SUFFIX => '.fastq' ); # temporary file with results pair1

#create log file
if ( -e $outputdir) {
	open (LOG, ">$outputdir/Log.txt") or die ("cannot create $outputdir/Log.txt");
}
else {
	open (LOG, ">Log.txt") or die ("cannot create Log.txt");
}

#determine the length of the first sequence before trimming
open (INPUT, $readspair1) or die "cannot open file $readspair1\n";
my $line = <INPUT>;
$line = <INPUT>;
chomp $line;
my $untrimmedlength = length $line; #length of sequences before trimming
close INPUT;

########### start cleanning ######################
print LOG datetime, " Initial count\n";
print LOG datetime, " File $readspair1, total FASTQ reads: ", count_fastq($readspair1), "\n"; # do a basic count
print LOG datetime, " File $readspair2, total FASTQ reads: ", count_fastq($readspair2), "\n"; # do a basic count
print LOG "\n";

#clip adapters
print LOG datetime, " Clipping adapters and quality filter with Trimmomatic\n";
clipadapters($readspair1, $readspair2);
print LOG datetime, " File with paired reads, FASTQ reads: ", count_fastq($paired_output), "\n"; 
print LOG datetime, " File with unpaired reads, FASTQ reads: ", count_fastq($unpaired_output), "\n";
print LOG "\n";

print LOG "\n";
#artifact filter
print LOG datetime, " Removing artifacts from unpaired file only\n";
filter_artifact($unpaired_output);
print LOG datetime, " File with unpaired reads, FASTQ reads: ", count_fastq($unpaired_output), "\n";
print LOG "\n";

#remove ribosomal sequence
print LOG datetime, " Removing ribosomal sequences using SortMeRNA\n";
`cp $paired_output t5`;
ribosome_removal($paired_output, 1);
print LOG datetime, " File with paired reads, FASTQ reads: ", count_fastq($paired_output), "\n";
ribosome_removal($unpaired_output, 0);
print LOG datetime, " File with unpaired reads, FASTQ reads: ", count_fastq($unpaired_output), "\n";

#print the data files
my $pairedoutname = $outputdir . "/" . $outputname .  "-paired.fq"; 
my $unpairedoutname = $outputdir . "/" . $outputname .  "-unpaired.fq";
`mv $paired_output $pairedoutname`; 
`mv $unpaired_output $unpairedoutname`;
print LOG datetime, " Paired and unpaired data file are written in $pairedoutname and $unpairedoutname\n";
print LOG "\n";

######## subroutines ###################
sub count_fastq {
#	print STDERR "count fastq...\n";
	my ($title) = @_;
	my $rawcount;

	my $txtcount = `wc -l $title`;
	if ($txtcount =~ /^(\d+)\s/) {
		my $count = $1/4;
		return (commify($count));
#		return(sprintf('%.3e', $count));
#		return ($count);
	}
	else {
		die "could not count lines in file $title using wc -l\n";
	}
}

sub clipadapters {
#	print STDERRR "clipping adapters...\n";	
	my ($inputfile1, $inputfile2) = @_;
	my $forward_paired = File::Temp->new( UNLINK => 1, SUFFIX => '.fastq' ); 
	my $reverse_paired = File::Temp->new( UNLINK => 1, SUFFIX => '.fastq' ); 
	my $forward_unpaired = File::Temp->new( UNLINK => 1, SUFFIX => '.fastq' );
	my $reverse_unpaired = File::Temp->new( UNLINK => 1, SUFFIX => '.fastq' ); 

	`java -jar $TRIMMOMATIC_PATH/trimmomatic-0.32.jar PE -threads $threads -phred33 $inputfile1 $inputfile2 $forward_paired $forward_unpaired $reverse_paired $reverse_unpaired ILLUMINACLIP:$TRIMMOMATIC_PATH/illuminaClipping.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:$MINLEN`;
	
	#merge the files
	open (INPUT1, $forward_paired) or die;
	open (INPUT2, $reverse_paired) or die;
	open (OUTPUT1, ">$paired_output") or die;
	while (my $line1 = <INPUT1>) {
		my $seq1 = $line1 . <INPUT1> . <INPUT1> . <INPUT1>;
		my $seq2 = <INPUT2> . <INPUT2> . <INPUT2> . <INPUT2>;
		print OUTPUT1 "$seq1";
		print OUTPUT1 "$seq2";
	}	
	`cat $forward_unpaired $reverse_unpaired > $unpaired_output`;
	close INPUT1;
	close INPUT2;
	close OUTPUT1; 
}

sub ribosome_removal{
	my ($infile, $paired) = @_;
	open (INPUT1, $infile) or die "cannot open file $infile\n";
        my $ribooutput = File::Temp->new( UNLINK => 1);
	my $ribooutput2 = File::Temp->new( UNLINK => 1);
	if ($paired) {
		`$SORTMELOC/sortmerna --ref $SORTMELOC/rRNA_databases/silva-arc-23s-id98.fasta,$SORTMELOC/index/silva-arc-23s-db:$SORTMELOC/rRNA_databases/misc_rRNA.fasta,$SORTMELOC/index/misc_rRNA-db:$SORTMELOC/rRNA_databases/silva-bac-16s-id90.fasta,$SORTMELOC/index/silva-bac-16s-db:$SORTMELOC/rRNA_databases/rfam-5.8s-database-id98.fasta,$SORTMELOC/index/rfam-5.8s-database-db:$SORTMELOC/rRNA_databases/silva-bac-23s-id98.fasta,$SORTMELOC/index/silva-bac-23s-db:$SORTMELOC/rRNA_databases/rfam-5s-database-id98.fasta,$SORTMELOC/index/rfam-5s-db:$SORTMELOC/rRNA_databases/silva-euk-18s-id95.fasta,$SORTMELOC/index/silva-euk-18s-db:$SORTMELOC/rRNA_databases/silva-euk-28s-id98.fasta,$SORTMELOC/index/silva-euk-28s-db:$SORTMELOC/rRNA_databases/silva-arc-16s-id95.fasta,$SORTMELOC/index/silva-arc-16s-db --reads $infile --other $ribooutput2 -a $threads --paired_in --fastx --sam --aligned $ribooutput`; 

		# this outputs the aligned ribosomes for Sandra
#		`cp $ribooutput aligned_paired-ribosomes`;

		#sortme still leaves unmatched pairs, this tries to remove them and put them into the unpaired file
		my $tempname = "$ribooutput2" . ".fastq"; #necessary because sortmeRNA adds .fastq to the output file without asking
		my $ref1;
		my $ref2;
		open (INPUT, $tempname) or die;
		open (OUTPUT1, ">$paired_output") or die;
		open (OUTPUT2, ">$unpaired_output") or die;
		(my $ref1) = <INPUT> =~ /\A([^\s]+)/; # get the first word, from http://stackoverflow.com/questions/4973229/perl-get-first-word-from-input-string
		my $ref2  = <INPUT> . <INPUT> . <INPUT>;
		while ((my $l1) = <INPUT> =~ /\A([^\s]+)/) {
			my $l2 = <INPUT> . <INPUT> . <INPUT>;
			if ($ref1 eq $l1) {
				print OUTPUT1 "$ref1\n";
				print OUTPUT1 $ref2;
				print OUTPUT1 "$l1\n";
				print OUTPUT1 $l2;
				$ref1 = $l1;
				$ref2 = $l2;
			}
			else {
				print OUTPUT2 "$ref1\n";
				print OUTPUT2 $ref2;
				$ref1 = $l1;
				$ref2 = $l2;
			}	
		}
		close INPUT;
		close OUTPUT1;
		close OUTPUT2;
	}
	else {
		`$SORTMELOC/sortmerna --ref $SORTMELOC/rRNA_databases/silva-arc-23s-id98.fasta,$SORTMELOC/index/silva-arc-23s-db:$SORTMELOC/rRNA_databases/misc_rRNA.fasta,$SORTMELOC/index/misc_rRNA-db:$SORTMELOC/rRNA_databases/silva-bac-16s-id90.fasta,$SORTMELOC/index/silva-bac-16s-db:$SORTMELOC/rRNA_databases/rfam-5.8s-database-id98.fasta,$SORTMELOC/index/rfam-5.8s-database-db:$SORTMELOC/rRNA_databases/silva-bac-23s-id98.fasta,$SORTMELOC/index/silva-bac-23s-db:$SORTMELOC/rRNA_databases/rfam-5s-database-id98.fasta,$SORTMELOC/index/rfam-5s-db:$SORTMELOC/rRNA_databases/silva-euk-18s-id95.fasta,$SORTMELOC/index/silva-euk-18s-db:$SORTMELOC/rRNA_databases/silva-euk-28s-id98.fasta,$SORTMELOC/index/silva-euk-28s-db:$SORTMELOC/rRNA_databases/silva-arc-16s-id95.fasta,$SORTMELOC/index/silva-arc-16s-db --reads $infile --other $ribooutput2 -a $threads --paired_in --fastx --sam --aligned $ribooutput`; 
		
		# this outputs the aligned ribosomes for Sandra
#               `cp $ribooutput aligned_unpaired-ribosomes`;

		my $tempname = "$ribooutput2" . ".fastq"; #necessary because sortmeRNA adds .fastq to the output file without asking
		`mv $tempname $unpaired_output`;
	}
}

sub filter_artifact {
	my ($inputfile) = @_;
	my $tx = File::Temp->new( UNLINK => 1, SUFFIX => '.fastq' ); # temporary file with hits
	`fastx_artifacts_filter -i $inputfile -o $tx -Q 33`;
	my $wctext = `wc -m $tx`; #word count of the ouptput used to check that fastx ran
	if ($wctext =~ /^0\s/) {
		die "There was an error, either fastx_artifacts_filter is not installed or all the sequences were removed as artifacts\n";
	}
	else {
		`mv $tx $inputfile`;
	}
}

# from perl cookbook, to put commas on big numbers
sub commify {
    my $text = reverse $_[0];
    $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $text
}

