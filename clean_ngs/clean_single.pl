#! /usr/bin/perl
# Wed 21 Mar 2012 02:30:34 PM PDT automates initial library prep
# Fri 23 Mar 2012 10:17:03 AM PDT modified to work with a pair of libraries
# Thu 19 Apr 2012 02:27:00 PM PDT updated to clean single file

# CONSTANTS
# adapter trimming parameters
my $TRIMMOMATIC_PATH = "./Trimmomatic-0.20";
my $RIBOSOME_BOWTIE2_FILE = "./arthropod_ribosomes";
my $MINLEN = 26;

#return date and time
sub datetime {
	use POSIX qw/strftime/;
	return (strftime('%D %T',localtime));
}

#### Main program
use strict;
require File::Temp;
use File::Temp ();
use File::Temp qw/ :seekable /;
use File::Basename;
use Getopt::Long;
use File::Path;

# other parameters
my $reads; #filename of first pair fastq
my $outputname; #base name for output files
my $outputdir; #outputdirectory
my $threads=1; #number of threads to use

#####read and check the inputs
GetOptions(
	'r:s'   => \$reads,
	'o:s'	=> \$outputname,
	'd:s'	=> \$outputdir,
	't:s'   => \$threads
);

unless ($reads) {
	die ("usage: perl clean_sinlge.pl -r <FASTQ file> -o <OPTIONAL output name> -d <OPTIONAL directory name for output> -t <OPTIONAL: number of threads to use, default $threads>\n");
}
unless ($outputname) {
	$outputname = basename($reads, ".fastq");
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
my $reads_output = File::Temp->new( UNLINK => 1, SUFFIX => '.fastq' ); # temporary file with chastity results

#create log file
if ( -e $outputdir) {
	open (LOG, ">$outputdir/Log.txt") or die ("cannot create $outputdir/Log.txt");
}
else {
	open (LOG, ">Log.txt") or die ("cannot create Log.txt");
}

#determine the length of the first sequence before trimming
open (INPUT, $reads) or die "cannot open file $reads\n";
my $line = <INPUT>;
$line = <INPUT>;
chomp $line;
my $untrimmedlength = length $line; #length of sequences before trimming
close INPUT;

########### start cleanning ######################
print LOG datetime, " Initial count\n";
print LOG datetime, " File $reads, total FASTQ reads: ", count_fastq($reads), "\n"; # do a basic count
print LOG "\n";

#clip adapters
print LOG datetime, " Clipping adapters and quality filter with Trimmomatic\n";
clipadapters($reads);
print LOG datetime, " File with reads, FASTQ reads: ", count_fastq($reads_output), "\n"; 
print LOG "\n";

#remove sequences that were not clipped
print LOG datetime, " Removing sequences that did NOT have a recognisable adapter sequence\n";
removeunclipped($reads_output, $untrimmedlength);
print LOG datetime, " After removing sequences that do not have adapter sequences left: ", count_fastq($reads_output), "\n";

# test to see if header of first line is compatible with chastity filter
open (my $fh, $reads_output) or die("ack -$!");
my $firstline = <$fh>;
if ($firstline =~ /:\S\s$/) {
	print LOG datetime, " Applying chastity filter\n";
	filter_chastity($reads_output);
	print LOG datetime, " File with reads, FASTQ reads: ", count_fastq($reads_output), "\n";
}
else {
	print LOG datetime, " Chastity filter not applied\n";
}
close $fh;

#
##chastity filter
#print LOG datetime, " Applying chastity filter\n";
#filter_chastity($reads_output);
#print LOG datetime, " File with reads, FASTQ reads: ", count_fastq($reads_output), "\n";
#print "\n";

#artifact filter
print LOG datetime, " Removing artifacts\n";
filter_artifact($reads_output);
print LOG datetime, " File with reads, FASTQ reads: ", count_fastq($reads_output), "\n";
print LOG "\n";

#remove ribosomal sequence
print LOG datetime, " Removing ribosomal sequences\n";
ribosome_removal($reads_output, 0);
print LOG datetime, " File with reads, FASTQ reads: ", count_fastq($reads_output), "\n";

#print the data files
my $readsoutname = $outputdir . "/" . $outputname .  ".fq"; 
`mv $reads_output $readsoutname`; 
print LOG datetime, " Processed file is written in $readsoutname\n";
print LOG "\n";

#print quality results
print LOG datetime, " Producing quality statistics\n";
print LOG datetime, " Quality boxplot and nucleotide distributions are in directory $outputdir in files: ", stats1($readsoutname), " \n";

sub clipadapters {
#       print STDERRR "clipping adapters...\n"; 
        my ($inputfile1) = @_;
        my $outputfile1 = File::Temp->new( UNLINK => 1, SUFFIX => '.fastq' );

	#run trimmomatic
        `java -classpath $TRIMMOMATIC_PATH/trimmomatic-0.20.jar org.usadellab.trimmomatic.TrimmomaticSE -threads $threads -phred33 $inputfile1 $reads_output ILLUMINACLIP:$TRIMMOMATIC_PATH/illuminaClipping.fa:2:40:15 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:$MINLEN`;
}

sub stats1 {
	my ($inputfile) = @_;
	my $basename = basename($inputfile, ".fq");
	my $ts = File::Temp->new( UNLINK => 1, SUFFIX => '.txt' ); # temporary file
	my $boxfilename = $basename . "-quality_boxplot.png";
	my $nucfilename = $basename . "-nuc_dist.png";
	`fastx_quality_stats -i $inputfile -o $ts -Q33`;
	`fastq_quality_boxplot_graph.sh -i $ts -o "$outputdir/$boxfilename" -t $basename`;
	`fastx_nucleotide_distribution_graph.sh -i $ts -o "$outputdir/$nucfilename" -t $basename`;
	return ($boxfilename . ", " . $nucfilename);
}

sub filter_chastity {
	my ($inputfile) = @_;
	open (INPUT, $inputfile) or die "cannot open file $inputfile\n";
	my $chasoutput = File::Temp->new( UNLINK => 1, SUFFIX => '.fastq' );
	open (OUTPUT, ">$chasoutput") or die "cannot open output file $chasoutput\n";
	while (my $line = <INPUT>) {
		if ($line =~ /^@/) {
			my $seq = $line . <INPUT> . <INPUT> . <INPUT>;
			if (($line =~ /:Y\s$/) || ($line =~ /:Y\/\d\s$/)) {
				print OUTPUT "$seq";
			}
		}
		else {
			die "file $inputfile does not match FASTQ at line:\n$line";
		}
	}
	close INPUT;
	close OUTPUT;
	`mv $chasoutput $reads_output`;
}

sub count_fastq {
#	print STDERR "count fastq...\n";
	my ($title) = @_;
	my $rawcount;

	my $txtcount = `wc -l $title`;
	if ($txtcount =~ /^(\d+)\s/) {
		my $count = $1/4;
		return(sprintf('%.3e', $count));
#		return ($count);
	}
	else {
		die "could not count lines in file $title using wc -l\n";
	}
}

sub ribosome_removal{
	my ($infile, $paired) = @_;
	open (INPUT1, $infile) or die "cannot open file $infile\n";
        my $ribooutput = File::Temp->new( UNLINK => 1, SUFFIX => '.fastq' );
	open (OUTPUT1, ">$ribooutput") or die;

	my $tx = File::Temp->new( UNLINK => 1, SUFFIX => '.sam' ); # temporary file with hits
	`bowtie2 -x $RIBOSOME_BOWTIE2_FILE -U $infile -S $tx -p $threads --local -k 1 --sam-nohead --sam-nosq`;

	#record the names of the files that matched
	open (INPUT2, $tx) or die "cannot open output of bowtie $tx";
	my %riboreads; #holds as key the reads that are ribosomes
	while (my $line = <INPUT2>) {
		if ($line =~ /^(\S+)\s\S+\s(\S+)\s/) {
			my $hit = $2;
			my $name = $1;
			unless ($hit =~ /\*/){
				$riboreads{$name} = 0;
			}
		}
		else {
			die "error reading sam line of file $tx at line\n$line";
		}
	}
	close INPUT2;

	#go through input file and report into output file only non-ribosomal sequences
	while (my $line1 = <INPUT1>) {
		my $seq1 = <INPUT1> . <INPUT1> . <INPUT1>;
		my $read2;
		if ($line1 =~ /^@(\S+)\s/) {
			$read2 = $1;
		}
		unless (exists $riboreads{$read2}) {
			print OUTPUT1 "$line1", "$seq1";
		}
	}
	
	`mv $ribooutput $reads_output`;

	close INPUT1;
	close OUTPUT1;
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

sub removeunclipped {
	my ($inputfile, $ulen) = @_;
	open (INPUT, $inputfile) or die "cannot open file $inputfile\n";
	my $outputfile = File::Temp->new( UNLINK => 1, SUFFIX => '.fastq' ); # temporary file with good reads
	open (OUTPUT, ">$outputfile") or die "cannot open output file $outputfile\n";

	while (my $l1 = <INPUT>) {
		my $l2 = <INPUT>;
		my $l3 = <INPUT>;
		my $l4 = <INPUT>;

		chomp $l2;
		if (length $l2 < $ulen) {
			print OUTPUT "$l1";
			print OUTPUT "$l2\n";
			print OUTPUT "$l3";
			print OUTPUT "$l4";
		}
	}

	close OUTPUT;
	`mv $outputfile $reads_output`; 
}
