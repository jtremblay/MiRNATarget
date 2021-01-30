#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
parseMiRNATargets.pl

PURPOSE:
From the 9C output of ssearch36, distribute miRNA into their site types. According to 
scheme found in PMC3499661

INPUT:
--infile <string>  : ssearch 9C out format file

OUTPUT:
STDOUT <string>    : standard output

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@nrc-cnrc.gc.ca

ENDHERE

## OPTIONS
my ($help, $infile);
my $verbose = 0;

GetOptions(
  'infile=s'                             => \$infile,
#  'outdir=s'                             => \$outdir,
  'verbose'                              => \$verbose,
  'help'                                 => \$help
);
if ($help) { print $usage; exit; }

die("--infile missing...\n") unless($infile);
#die("--outdir_fasta missing...\n") unless($outdir);

## MAIN
my %hash;
open(IN, "<".$infile) or die "Can't open $infile\n";
print STDERR "Processing $infile\n";
my $curr_query = "";
my $curr_target = "";
my $curr_query_length;
my $counter_match_nucl_string = 0;
my $curr_start;
my $curr_end;
while(<IN>){
    chomp;
    if($_ =~ m/^#/){
        next;
    }
    if($_ =~ m/^>>><<</){
        $counter_match_nucl_string = 0;
        $curr_start = 0;
        $curr_end = 0;
        next;
    }
    if($_ =~ m/^>>>(\S+), (\d+) nt/){
        $curr_query = $1;
        $curr_query_length = $2;
        $hash{$curr_query}{length} = $curr_query_length;
        next;
    }
    if($_ =~ m/^>>(\S+) .*\(\d+ nt\)$/){
        $curr_target = $1;
        #$curr_target =~ s/\s+$//;
        #print STDERR $curr_target."\n";
        $hash{$curr_query}{$curr_target}{name} = $curr_target;
        next;
    }
    #GL635794.1 Neisseria mucosa C102 genomic scaffold superc (494621) [f]
    if($_ =~ m/^\S+\s.*\(\d+\) \[/){
        #print STDERR $_."\n";
        next;
    }
    if($_ =~ m/^\S+\s+([ACGTU-]*)\s*$/ && $counter_match_nucl_string == 0){
        # Found a query nucl alignment string.
        $curr_start = $-[1];
        $curr_end = $+[1];
        #print STDERR "curr_start: ".$curr_start."\n";
        #print STDERR "curr_end: ".$curr_end."\n";
        my $str = substr($_, $curr_start, $curr_end - $curr_start);
        $hash{$curr_query}{$curr_target}{query_aln} = $str;
        
        $counter_match_nucl_string = 1;
        next;
    }
    if($_ =~ m/^\s+[\.\:]/ && $counter_match_nucl_string == 1){
        # Found a match string.
        #extract substring at previously found positions.
        my $str = substr($_, $curr_start, $curr_end - $curr_start);
        $hash{$curr_query}{$curr_target}{match_aln} = $str;
        next;
    }
    #Smith-Waterman score: 279; 89.5% identity (100.0% similar) in 19 nt overlap (2-20:52164-52182)
    if($_ =~ m/^Smith-Waterman.*\(\d+-\d+:(\d+)-(\d+)\)$/){
        my $start = $1;
        my $end = $2;
        #print STDERR "1: ".$start."\n";
        #print STDERR "2: ".$end."\n";
        $hash{$curr_query}{$curr_target}{start} = $start;
        $hash{$curr_query}{$curr_target}{end} = $end;
        next;
    }
    if($_ =~ m/^\S+\s+[ACGTUYNWRMKS-]*\s*$/ && $counter_match_nucl_string == 1){
        #extract substring at previously found positions.
        my $str = substr($_, $curr_start, $curr_end - $curr_start);
        # Found subject string of current match.
        $hash{$curr_query}{$curr_target}{subject_aln} = $str;
       
        if(!exists($hash{$curr_query}{$curr_target}{query_aln})){
            print STDERR "line # ".$.."\n";
            print STDERR Dumper(\%hash);
            die "Problem at $curr_query\n$curr_target\n";
        }

        $counter_match_nucl_string = 0;
        next;
    }
}
close(IN);

print STDERR Dumper(\%hash);
print STDERR "Done parsing ssearch 9C outfmt file. Will now parse each alignments...\n";

# Finally sort and print hash.
my %hash_aln;
foreach my $query (keys %hash) {
    #my $length = $hash{$query}{length};

    foreach my $subject (keys %{ $hash{$query} }) {
        next if($subject eq "length");
        my $aln_str = reverse($hash{$query}{$subject}{match_aln});
        my @aln_char = split(//, $aln_str);
        my $length = scalar(@aln_char);

        $hash_aln{reverse($aln_str)}++;
        
        #my $i = 0;
        #my $number_of_matches = 0;
        #for my $el (@aln_char){
        #    if($el eq ":"){
        #        $number_of_matches++;
        #    }
        #}
    }
}

# Take 20 most abundant matches.
my %hash_aln_ma;
#my $i = 0;
foreach my $match_string (sort { $hash_aln{$b} <=> $hash_aln{$a} } keys %hash_aln) {
    #next if($i > 200);
    $hash_aln_ma{$match_string}{freq} = $hash_aln{$match_string};
    print STDERR $match_string."=>".$hash_aln{$match_string}."\n";
    #$i++;
}
print STDERR Dumper(\%hash_aln_ma);

foreach my $query (keys %hash) {

    foreach my $subject (keys %{ $hash{$query} }) {
        next if($subject eq "length");
        my $aln_str = $hash{$query}{$subject}{match_aln};
        if(exists $hash_aln_ma{$aln_str}){
            print STDOUT $query."\t".$subject."\t".$aln_str."\t".$hash{$query}{$subject}{query_aln}."\t".$hash{$query}{$subject}{subject_aln}."\t".$hash{$query}{$subject}{start}."\t".$hash{$query}{$subject}{end}."\n";
        }
    }
}


# 

exit;


#close(OUT_FULLMATCH);
#close(OUT_P1);
#close(OUT_P2);
#close(OUT_P3);
#close(OUT_P4);
#close(OUT_P5);
#open(OUT_FULLMATCH, ">".$outdir."/fullmatch.txt") or die "Can't open ".$outdir."/fullmatch.txt";
#open(OUT_P1, ">".$outdir."/P1.txt") or die "Can't open ".$outdir."/P1.txt";
#open(OUT_P2, ">".$outdir."/P2.txt") or die "Can't open ".$outdir."/P2.txt";
#open(OUT_P3, ">".$outdir."/P3.txt") or die "Can't open ".$outdir."/P3.txt";
#open(OUT_P4, ">".$outdir."/P4.txt") or die "Can't open ".$outdir."/P4.txt";
#open(OUT_P5, ">".$outdir."/P5.txt") or die "Can't open ".$outdir."/P5.txt";
#        if($number_of_matches == $length){
#            print OUT_FULLMATCH ">>".$query."\n";
#            print OUT_FULLMATCH ">".$subject."\n";
#            print OUT_FULLMATCH $hash{$query}{$subject}{query_aln}."\n";
#            print OUT_FULLMATCH $hash{$query}{$subject}{match_aln}."\n";
#            print OUT_FULLMATCH $hash{$query}{$subject}{subject_aln}."\n\n";
#        }else{
#            #case 1
#            if( 
#                $aln_char[0]  eq " " &&
#                $aln_char[1]  eq ":" &&
#                $aln_char[2]  eq ":" &&
#                $aln_char[3]  eq ":" &&
#                $aln_char[4]  eq ":" &&
#                $aln_char[5]  eq ":" &&
#                $aln_char[6]  eq ":" &&
#                $aln_char[7]  eq ":" &&
#                $aln_char[8]  eq " " &&
#                $aln_char[9]  eq " " &&
#                $aln_char[10] eq " " &&
#                $aln_char[11] eq " " &&
#                $aln_char[12] eq " " &&
#                $aln_char[13] eq " " &&
#                $aln_char[14] eq " " &&
#                $aln_char[15] eq " " &&
#                $aln_char[16] eq " " &&
#                $aln_char[17] eq "."){
#                
#                print OUT_P1 ">>".$query."\n";
#                print OUT_P1 ">".$subject."\n";
#                print OUT_P1 $hash{$query}{$subject}{query_aln}."\n";
#                print OUT_P1 $hash{$query}{$subject}{match_aln}."\n";
#                print OUT_P1 $hash{$query}{$subject}{subject_aln}."\n\n";
#            }elsif( 
#                $aln_char[0]  eq " " &&
#                $aln_char[1]  eq ":" &&
#                $aln_char[2]  eq ":" &&
#                $aln_char[3]  eq ":" &&
#                $aln_char[4]  eq ":" &&
#                $aln_char[5]  eq "." &&
#                $aln_char[6]  eq ":" &&
#                $aln_char[7]  eq ": " &&
#                $aln_char[8]  eq " " &&
#                $aln_char[9]  eq ":" &&
#                $aln_char[10] eq " " &&
#                $aln_char[11] eq ":" &&
#                $aln_char[12] eq ":" &&
#                $aln_char[13] eq ":" &&
#                $aln_char[14] eq ":" &&
#                $aln_char[15] eq ":" &&
#                $aln_char[16] eq ":" &&
#                $aln_char[17] eq ":" &&
#                $aln_char[18] eq ":" &&
#                $aln_char[19] eq ":" &&
#                $aln_char[20] eq "."){
#                
#                print OUT_P2 ">>".$query."\n";
#                print OUT_P2">".$subject."\n";
#                print OUT_P2 $hash{$query}{$subject}{query_aln}."\n";
#                print OUT_P2 $hash{$query}{$subject}{match_aln}."\n";
#                print OUT_P2 $hash{$query}{$subject}{subject_aln}."\n\n";
#            }elsif( 
#                $aln_char[0]  eq ":" &&
#                $aln_char[1]  eq ":" &&
#                $aln_char[2]  eq ":" &&
#                $aln_char[3]  eq ":" &&
#                $aln_char[4]  eq "." &&
#                $aln_char[5]  eq ":" &&
#                $aln_char[6]  eq ":" &&
#                $aln_char[7]  eq ":" &&
#                $aln_char[8]  eq ":" &&
#                $aln_char[9]  eq ":" &&
#                $aln_char[10] eq ":" &&
#                $aln_char[11] eq ":" &&
#                $aln_char[12] eq ":" &&
#                $aln_char[13] eq ":" &&
#                $aln_char[14] eq ":" &&
#                $aln_char[15] eq ":" &&
#                $aln_char[16] eq ":" &&
#                $aln_char[17] eq ":" &&
#                $aln_char[18] eq ":" &&
#                $aln_char[19] eq ":" &&
#                $aln_char[20] eq ":"){
#                
#                print OUT_P3 ">>".$query."\n";
#                print OUT_P3 ">".$subject."\n";
#                print OUT_P3 $hash{$query}{$subject}{query_aln}."\n";
#                print OUT_P3 $hash{$query}{$subject}{match_aln}."\n";
#                print OUT_P3 $hash{$query}{$subject}{subject_aln}."\n\n";
#            }elsif( 
#                $aln_char[0]  eq " " &&
#                $aln_char[1]  eq " " &&
#                $aln_char[2]  eq "." &&
#                $aln_char[3]  eq " " &&
#                $aln_char[4]  eq ":" &&
#                $aln_char[5]  eq ":" &&
#                $aln_char[6]  eq ":" &&
#                $aln_char[7]  eq ":" &&
#                $aln_char[8]  eq ":" &&
#                $aln_char[9]  eq ":" &&
#                $aln_char[10] eq ":" &&
#                $aln_char[11] eq ":" &&
#                $aln_char[12] eq ":" &&
#                $aln_char[13] eq ":" &&
#                $aln_char[14] eq ":" &&
#                $aln_char[15] eq ":" &&
#                $aln_char[16] eq " " &&
#                $aln_char[17] eq "." &&
#                $aln_char[18] eq " " &&
#                $aln_char[19] eq " " &&
#                $aln_char[20] eq " "){
#                
#                print OUT_P4 ">>".$query."\n";
#                print OUT_P4 ">".$subject."\n";
#                print OUT_P4 $hash{$query}{$subject}{query_aln}."\n";
#                print OUT_P4 $hash{$query}{$subject}{match_aln}."\n";
#                print OUT_P4 $hash{$query}{$subject}{subject_aln}."\n\n";
#            }elsif( 
#                $aln_char[0]  eq "." &&
#                $aln_char[1]  eq ":" &&
#                $aln_char[2]  eq ":" &&
#                $aln_char[3]  eq ":" &&
#                $aln_char[4]  eq ":" &&
#                $aln_char[5]  eq " " &&
#                $aln_char[6]  eq ":" &&
#                $aln_char[7]  eq ":" &&
#                $aln_char[8]  eq ":" &&
#                $aln_char[9]  eq " " &&
#                $aln_char[10] eq " " &&
#                $aln_char[11] eq " " &&
#                $aln_char[12] eq "." &&
#                $aln_char[13] eq ":" &&
#                $aln_char[14] eq " " &&
#                $aln_char[15] eq " " &&
#                $aln_char[16] eq ":" &&
#                $aln_char[17] eq " " &&
#                $aln_char[18] eq ":" &&
#                $aln_char[19] eq " " &&
#                $aln_char[20] eq " "){
#                
#                print OUT_P5 ">>".$query."\n";
#                print OUT_P5 ">".$subject."\n";
#                print OUT_P5 $hash{$query}{$subject}{query_aln}."\n";
#                print OUT_P5 $hash{$query}{$subject}{match_aln}."\n";
#                print OUT_P5 $hash{$query}{$subject}{subject_aln}."\n\n";
#            }
#        }
