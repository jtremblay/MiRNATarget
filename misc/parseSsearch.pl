#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;
use File::Basename;

my $usage=<<'ENDHERE';
NAME:
parseSsearch.pl

PURPOSE:
From the default output of ssearch36, this script parses each alignment to a .tsv format that can then be used in the companion 
script parseMiRNATargets.pl.

INPUT:
--infile <string>                : ssearch default out format file

OUTPUT:
STDOUT <string>                  : standard output. Alignments in tsv format.

NOTES:
It is not entirely clear which ssearch36 parameters are used by psRNATarget, but using the following SSEARCH followed by the execution of this script with default parameters gives identical results to psRNATarget.
ssearch36 -f -8 -g -3 -E 10000 -T 8 -b 200 -r +4/-3 -n -U -W 10 -N 20000 -i <input_mirna.fasta> <reference_genome.fasta> > <output_file>
If you provide an already rev-complemented fasta file, you can omit the -i argument.

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@nrc-cnrc.gc.ca

ENDHERE

## OPTIONS
my ($help, $infile, $rev);
my $verbose = 0;

GetOptions(
  'infile=s'                   => \$infile,
  'rev'                        => \$rev,
  'verbose'                    => \$verbose,
  'help'                       => \$help
);
if ($help) { print $usage; exit; }
print STDERR "##########################################\n";
print STDERR "# Running parseSsearch.pl ...            #\n";
print STDERR "##########################################\n";

# open infile and start parsing alignments.
open(IN, "<", $infile) or die "Can't open $infile\n";
print STDERR "Processing $infile\n";
my %hash;
my $curr_query = "";
my $curr_target = "";
my $curr_aln = "";
my $curr_query_length;
my $counter_match_nucl_string = 0;
my $curr_start;
my $curr_end;
my $start;
my $end;
my $hsp;
my $q_end;
my $q_start;
my $query_str;
my $aln_str;
my $subject_str;
my $i = 0;
my $strand;
while(<IN>){
    chomp;
    if($_ =~ m/^#/){
        next;
    }
    if($_ =~ m/\d+>>>(\S+) - (\d+) nt/){
        if($verbose){ print STDERR "\n---------------------------------------------\n"; }
        $counter_match_nucl_string = 0;
        $curr_start = 0;
        $curr_end = 0;
        $curr_query = $1;
        $curr_query_length = $2;
        $i = 0;

        if($verbose){
            print STDERR "curr_query: ".$curr_query."\n";
        }

        next;
    }
    if($_ =~ m/^>>(\S+) .*\((\d+) nt\)$/){
        $i++;
        $curr_target = $1."_".$i;
        next;
    }
    #Smith-Waterman score: 279; 89.5% identity (100.0% similar) in 19 nt overlap (2-20:52164-52182)
    if($_ =~ m/^Smith-Waterman.*\((\d+)-(\d+):(\d+)-(\d+)\)$/){
        my $diff_start = $1 - 1;
        my $diff_end = $curr_query_length;
        $start = $3 - $diff_start;
        $end = $4 + $diff_end;
        $hsp = $4 - $3;
        if($rev){
            $q_start = $1;
            $q_end = $2;
            $strand = "-";
        }else{
            # If fwd, qstart is in reverse orientation.
            $q_start = $2;
            $q_end = $1;
            $strand = "+";
        }
        next;
    }
    if($_ =~ m/^>--$/){ #same contig, but in another location
        $i++;
        next;
    }
    if($_ =~ m/^\S+\s+([ACGTU-]*)\s*$/ && $counter_match_nucl_string == 0){
        # Found a query nucl alignment string.
        $curr_start = $-[1];
        $curr_end = $+[1];
        my $str = substr($_, $curr_start, $curr_end - $curr_start);
        $hash{$curr_query}{$curr_target."_".$i}{query_aln} = 1;
        $query_str = $str;
        
        $counter_match_nucl_string = 1;
        next;
    }
    if($_ =~ m/^\s+[.:]/ && $counter_match_nucl_string == 1){
        # Found a match string.
        #extract substring at previously found positions.
        my $str = substr($_, $curr_start, $curr_end - $curr_start);
        $aln_str = $str;
        next;
    }
    if($_ =~ m/^\S+\s+[ACGTBDHUYNVWRMKS-]*\s*$/ && $counter_match_nucl_string == 1){
        #extract substring at previously found positions.
        my $str = substr($_, $curr_start, $curr_end - $curr_start);
        # Found subject string of current match.
        $subject_str = $str;
       
        if(!exists($hash{$curr_query}{$curr_target."_".$i}{query_aln})){
            die "Problem at $curr_query\n$curr_target"."_".$i."\n".$_."\n"."line number ".$.." in the file\n";
        }

        $counter_match_nucl_string = 0;

        print STDOUT $curr_query."\t".$curr_target."\t".$aln_str."\t".$query_str."\t".$subject_str."\t".$q_start."\t".$q_end."\t".$start."\t".$end."\t".$strand."\t".$hsp."\n";
        next;
    }
}
close(IN);
print STDERR  Dumper(\%hash) if($verbose);

print STDERR "Done parsing ssearch standard format outfmt file.\n";

sub complement_IUPAC {
   my $dna = shift;

   # complement the reversed DNA sequence
   my $comp = $dna;
   $comp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy-/TVGHCDKNYSAABWXRtvghcdknysaabwxr-/;
   return $comp;
}

exit;


