#!/usr/bin/perl
# Wed 21 Mar 2012 02:30:34 PM PDT automates initial library prep
# Fri 23 Mar 2012 10:17:03 AM PDT modified to work with a pair of libraries

use strict;
require File::Temp;
use File::Temp ();
use File::Temp qw/ :seekable /;
use File::Basename;
use Getopt::Long;
use File::Path;
use common_scripts;

# adapter trimming parameters
my $TRIMMOMATIC_PATH = "./Trimmomatic-0.20";
my $RIBOSOME_BOWTIE2_FILE = "/home/peter/Work/ribosomal_database/arthropod_ribosomes";
my $MINLEN = 20;

# other parameters
my $readspair1; #filename of first pair fastq
my $readspair2; #filename of second pair fastq file
my $outputname; #base name for output files
my $outputdir; #outputdirectory

#####read and check the inputs
GetOptions(
	'1:s'   => \$readspair1,
	'2:s'	=> \$readspair2,
	'o:s'	=> \$outputname,
	'd:s'	=> \$outputdir,
);

unless ($readspair1 and $readspair2) {
	die ("usage: perl clean_pair.pl -1 <FASTQ file pair 1> -2 <FASTQ file pair 2> -o <OPTIONAL output name> -d <OPIONAL output directory (default current directory)>\n");
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
my $paired_output = File::Temp->new( UNLINK => 1, SUFFIX => '.fastq' ); # temporary file with chastity results pair1
my $unpaired_output = File::Temp->new( UNLINK => 1, SUFFIX => '.fastq' ); # temporary file with chastity results pair1


########### start cleanning ######################
print datetime, " Initial count\n";
print datetime, " File $readspair1, total FASTQ reads: ", count_fastq($readspair1), "\n"; # do a basic count
print datetime, " File $readspair2, total FASTQ reads: ", count_fastq($readspair2), "\n"; # do a basic count
print "\n";

#clip adapters
print datetime, " Clipping adapters and quality filter with Trimmomatic\n";
clipadapters($readspair1, $readspair2);
print datetime, " File with paired reads, FASTQ reads: ", count_fastq($paired_output), "\n"; 
print datetime, " File with unpaired reads, FASTQ reads: ", count_fastq($unpaired_output), "\n";
print "\n";

#chastity filter
print datetime, " Applying chastity filter\n";
filter_chastity_pair($paired_output);
print datetime, " File with paired reads, FASTQ reads: ", count_fastq($paired_output), "\n";
filter_chastity($unpaired_output);
print datetime, " File with unpaired reads, FASTQ reads: ", count_fastq($unpaired_output), "\n";
print "\n";

#artifact filter
print datetime, " Removing artifacts from unpaired file only\n";
filter_artifact($unpaired_output);
print datetime, " File with unpaired reads, FASTQ reads: ", count_fastq($unpaired_output), "\n";
print "\n";

#remove ribosomal sequence
print datetime, " Removing ribosomal sequences\n";
ribosome_removal($paired_output, 1);
print datetime, " File with paired reads, FASTQ reads: ", count_fastq($paired_output), "\n";
ribosome_removal($unpaired_output, 0);
print datetime, " File with unpaired reads, FASTQ reads: ", count_fastq($unpaired_output), "\n";

#print the data files
my $pairedoutname = $outputdir . "/" . $outputname .  "-paired.fq"; 
my $unpairedoutname = $outputdir . "/" . $outputname .  "-unpaired.fq";
`mv $paired_output $pairedoutname`; 
`mv $unpaired_output $unpairedoutname`;
print datetime, " Paired and unpaired data file are written in $pairedoutname and $unpairedoutname\n";
print "\n";

#print quality results
print datetime, " Producing quality statistics\n";
print datetime, " Paired file quality boxplot and nucleotide distributions are in directory $outputdir in files: ", stats1($pairedoutname), " \n";
print datetime, " Unaired file quality boxplot and nucleotide distributions are in directory $outputdir in files: ", stats1($unpairedoutname), " \n";


######## subroutines ###################

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

sub filter_chastity_pair {
	my ($infile) = @_;
	open (INPUT, $infile) or die "cannot open file $infile\n";
        my $chasoutput = File::Temp->new( UNLINK => 1, SUFFIX => '.fastq' );
	open (OUTPUT, ">$chasoutput") or die;
	open (OUTPUT2, ">>$unpaired_output") or die;

	while (my $line = <INPUT>) { #not testing the file here, assuming it's fastq
		my $seq1 = $line. <INPUT> . <INPUT> . <INPUT>; #forward sequence
		my $seq2 = <INPUT>. <INPUT> . <INPUT> . <INPUT>; #reverse sequence
		
		if (($seq1 =~ /^\S+:Y\s/) and ($seq2 =~ /^\S+:N\s/) ) {
			print OUTPUT2 "$seq1";
		}
		elsif (($seq1 =~ /^\S+:N\s/) and ($seq2 =~ /^\S+:Y\s/) ) {
			print OUTPUT2 "$seq2";
		}
		elsif (($seq1 =~ /^\S+:Y\s/) and ($seq2 =~ /^\S+:Y\s/) ) {
			print OUTPUT "$seq1";
			print OUTPUT "$seq2";
		}
		elsif (($seq1 =~ /^\S+:N\s/) and ($seq2 =~ /^\S+:N\s/) ) {
		
		}
		else {
			die "error reading files for chastity filter\n$seq1\n\n$seq2";
		}
	}
	close INPUT;
	close OUTPUT;
	close OUTPUT2;
	`mv $chasoutput $paired_output`;
}

sub filter_chastity {
	my ($inputfile) = @_;
	open (INPUT, $inputfile) or die "cannot open file $inputfile\n";
	my $chasoutput = File::Temp->new( UNLINK => 1, SUFFIX => '.fastq' );
	open (OUTPUT, ">$chasoutput") or die "cannot open output file $chasoutput\n";
	while (my $line = <INPUT>) {
		if ($line =~ /^@/) {
			my $seq = $line . <INPUT> . <INPUT> . <INPUT>;
			if ($line =~ /:Y\s$/) {
				print OUTPUT "$seq";
			}
		}
		else {
			die "file $inputfile does not match FASTQ at line:\n$line";
		}
	}
	close INPUT;
	close OUTPUT;
	`mv $chasoutput $unpaired_output`;
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

sub clipadapters {
#	print STDERRR "clipping adapters...\n";	
	my ($inputfile1, $inputfile2) = @_;
	my $forward_paired = File::Temp->new( UNLINK => 1, SUFFIX => '.fastq' ); 
	my $reverse_paired = File::Temp->new( UNLINK => 1, SUFFIX => '.fastq' ); 
	my $forward_unpaired = File::Temp->new( UNLINK => 1, SUFFIX => '.fastq' );
	my $reverse_unpaired = File::Temp->new( UNLINK => 1, SUFFIX => '.fastq' ); 

	`java -classpath $TRIMMOMATIC_PATH/trimmomatic-0.20.jar org.usadellab.trimmomatic.TrimmomaticPE -phred33 $inputfile1 $inputfile2 $forward_paired $forward_unpaired $reverse_paired $reverse_unpaired ILLUMINACLIP:$TRIMMOMATIC_PATH/illuminaClipping.fa:2:40:15 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:$MINLEN`;

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
        my $ribooutput = File::Temp->new( UNLINK => 1, SUFFIX => '.fastq' );
	open (OUTPUT1, ">$ribooutput") or die;

	my $tx = File::Temp->new( UNLINK => 1, SUFFIX => '.sam' ); # temporary file with hits
	`bowtie2 -x $RIBOSOME_BOWTIE2_FILE -U $infile -S $tx -p 8 --local -k 1 --sam-nohead --sam-nosq`;

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
	if ($paired) { #paired case
		while (my $line1 = <INPUT1>) {
			my $seq1 = <INPUT1> . <INPUT1> . <INPUT1>;
			my $seq2 = <INPUT1> . <INPUT1> . <INPUT1> . <INPUT1>;
	
			my $read1 = substr $line1, 1;
			chomp $read1;
	
			unless (exists $riboreads{$read1}) {
				print OUTPUT1 "$line1", "$seq1", "$seq2";
			}	
		}
		`mv $ribooutput $paired_output`; 
	}
	else { #unpaired
		while (my $line1 = <INPUT1>) {
			my $seq1 = <INPUT1> . <INPUT1> . <INPUT1>;
			my $read1 = substr $line1, 1;
			chomp $read1;
			unless (exists $riboreads{$read1}) {
				print OUTPUT1 "$line1", "$seq1";
			}
		}
		
		`mv $ribooutput $unpaired_output`;
	}
	close INPUT1;
	close OUTPUT1;
}

sub filter_artifact {
	my ($inputfile) = @_;
	my $tx = File::Temp->new( UNLINK => 1, SUFFIX => '.fastq' ); # temporary file with hits
	`fastx_artifacts_filter -i $inputfile -o $tx -Q 33`;
	`mv $tx $inputfile`;
}
