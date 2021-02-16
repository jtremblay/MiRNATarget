#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
mergeMiRNATargets.pl

PURPOSE:
merge to output of parseMiRNATargets.pl both fwd and rev into a 
single file.

INPUT:
--infile_fwd <string>  : infile fwd hits
--infile_rev <string>  : infile rev hits

OUTPUT:
STDOUT <string>        : standard output

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@nrc-cnrc.gc.ca

ENDHERE

## OPTIONS
my ($help, $infile_fwd, $infile_rev);
my $verbose = 0;

GetOptions(
  'infile_fwd=s' => \$infile_fwd,
  'infile_rev=s' => \$infile_rev,
  'verbose'      => \$verbose,
  'help'         => \$help
);
if ($help) { print $usage; exit; }

die("--infile_fwd missing...\n") unless($infile_fwd);
die("--infile_rev missing...\n") unless($infile_rev);

## MAIN
my %hash;
open(IN, "<".$infile_fwd) or die "Can't open $infile_fwd\n";
while(<IN>){
    chomp;
    my @row_orig = split(/\t/, $_);
    my @row = split(/\t/, $_);
    my $query = shift(@row);
    my $subject = shift(@row);
    $hash{$query}{$subject}{fwd} = join("\t", @row_orig);
}
close(IN);

open(IN, "<".$infile_rev) or die "Can't open $infile_rev\n";
while(<IN>){
    chomp;
    my @row_orig = split(/\t/, $_);
    my @row = split(/\t/, $_);
    my $query = shift(@row);
    my $subject = shift(@row);
    #my $aln_str = shift(@row);
    #my $query_str = shift(@row);
    #my $subject_str = shift(@row);
    #my $end = shift(@row);
    #my $start = shift(@row);
    #my $evalue = shift(@row);

    #$hash{$query}{$subject}{rev} = "$query\t$subject\t$aln_str\t$query_str\t$subject_str\t$start\t$end\t$evalue\t-";
    $hash{$query}{$subject}{rev} = join("\t", @row_orig);
}
close(IN);

foreach my $query (keys %hash) {
    foreach my $subject (keys %{ $hash{$query} }) {
        foreach my $strand (keys %{ $hash{$query}{$subject} }) {
            my $line = $hash{$query}{$subject}->{$strand};
            print STDOUT $line."\n";
        }
    }
}
