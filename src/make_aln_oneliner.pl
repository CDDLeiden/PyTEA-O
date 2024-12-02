#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;
use File::Path qw(make_path);

my $name = "make_aln_oneliner.pl";
my $version = "0.0.1";
my $updated = "06-08-2024";

my $usage = <<EXIT;

NAME		${name}
VERSION		${version}
UPDATED		${updated}

SYNOPSIS	Transforms windowed alignment files created with
		clustal to one-liner alignments for reading with TEA-O

EXIT

my @infiles;
my $outdir = "./one_liners";
GetOptions(
	"i|input=s{1,}" => \@infiles,
	"o|outdir=s" => \$outdir,
);

unless (-d $outdir){
	make_path($outdir,{mode=>0755});
}


foreach my $infile (@infiles){

	my ($filename,$dir,$ext) = fileparse($infile,'\..+');
	my $tempdir = $outdir."/".$filename;
	my $outfile = $tempdir."/".$filename.".olaln";

	unless (-d $tempdir){
		make_path($tempdir,{mode=>0755}) or die("Unable to create directory $tempdir: $!\n")
	}

	open IN, "<", $infile or die("Unable to open alignment file $infile: $!\n\n");

	my @loci;
	my %alignments;
	my @sequence;
	while(my $line = <IN>){
		chomp($line);
		
		if ($line =~ /^(\S+(?:\.\S+)?)\s+(\S+)/){

			if ($1 eq 'CLUSTAL'){
				next;
			}

			if (!$alignments{$1}){
				push(@loci,$1);
			}
			$alignments{$1} .= $2;
		}
	}

	close IN;

	open OUT, ">", $outfile or die("Unable to open output file $outfile: $!\n\n");

	foreach my $locus (@loci){

		my $alignment = $alignments{$locus};

		print OUT ">".$locus."\n".$alignment."\n";

		(my $seq = $alignment) =~ s/-//g;

		@sequence = unpack("(A60)*",$seq);

	}

	close OUT;

	my $fastadir = $tempdir."/FASTA";

	unless (-d $fastadir){
		make_path($fastadir,{mode=>0755}) or die("Unable to create directory $fastadir: $!\n")
	}

	foreach my $locus (@loci){

		my $outfile = $fastadir."/".$locus.".fasta";

		open OUT, ">", $outfile or die("Unable to open output file $outfile: $!\n\n");

		print OUT ">".$locus."\n";

		foreach my $line (@sequence){
			print OUT $line."\n";
		}


		close OUT;
	}


}
