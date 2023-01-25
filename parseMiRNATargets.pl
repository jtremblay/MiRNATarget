#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
parseMiRNATargets.pl

PURPOSE:
From the default output of ssearch36, distribute miRNA into their site types. According to the 
scheme found in PMC3499661. Results of this script matches almost exactly the ones outputted 
by psRNATarget. The scoring sheme gives exactly the same results as psRNATarget.

INPUT:
--infile <string>                : ssearch default out format file
--E_cutoff <float>               : default = 5.0  - Expectation value is the penalty for the mismatches 
                                   between mature small RNA and the target sequence. A higher value indicates 
                                   less similarity (and possibility) between small RNA and the target candidate. 
                                   The default penalty rule is set up by the scoring schema. Maximum expectation is the cutoff; 
                                   any small RNA-target pair with an expectation less than the cutoff will be discarded in the final result. 
                                   The recommended values are 3.0-5.0 depending on the scoring schema. 
--penalty_multiplier <string>    : default = 1.5. In the seed region (2-13 nt from the 5' of the miRNA strand), multiply mismatches by 
                                   penalty_multiplier. Only mismatches are multiplied, not G:U pairs and/or the '.' alignment
                                   caracters given by SSEARCH.
--num_mismatch_seed <int>        : default = 2 - Maximum of allowed mismatches in the seed region, excluding G:U pairs.
--hsp_cutoff <int>               : default = 17 - HSPs shorter than this value will be discarded.
--gap_cutoff <int>               : default = 1 - alignments having more than <int> gaps willbe discarded.
--total_mismatches_cutoff <int>  : default = 8 - alignments showing more than <int> mismatches will be discarded.
--GUs_cutoff <int>               : default = 6 - alignments having more than <int> mismatches will be discarded.
--keep_target_suffix             : set flag if you wish to keep the temporary _<int> suffix appended at the end
                                   of subjects IDs.
--rev                            : Set flag if alignments of miRNAs were done on revcomp subject sequences.
--verbose                        : Set flag for debugging.
--maximum_alignment_length <int> : default = 21 - alignments longer than this value will be discarded.

OUTPUT:
STDOUT <string>    : standard output

NOTES:
It is not entirely clear which ssearch36 parameters are used by psRNATarget, but using the following parameters gives identical results to psRNATarget.
ssearch36 -f -16 -g -10.0 -E 200 -T 8 -b 200 -r +4/-3 -n -U -W 10 -N 1000 -i <input_mirna.fasta> <reference_genome.fasta> > <output_file>
If you provide an already rev-complemented fasta file, you can omit the -i argument.

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@nrc-cnrc.gc.ca

ENDHERE

## OPTIONS
my ($help, $infile, $E_cutoff, $num_mismatch_seed, $hsp_cutoff, $gap_cutoff, $total_mismatches_cutoff, $GUs_cutoff, $keep_target_suffix, $penalty_multiplier, $rev, $maximum_alignment_length);
my $verbose = 0;

GetOptions(
  'infile=s'                   => \$infile,
  'E_cutoff=f'                 => \$E_cutoff,
  'num_mismatch_seed=i'        => \$num_mismatch_seed,
  'hsp_cutoff=i'               => \$hsp_cutoff,
  'gap_cutoff=i'               => \$gap_cutoff,
  'GUs_cutoff=i'               => \$GUs_cutoff,
  'penalty_multiplier=f'       => \$penalty_multiplier,
  'rev'                        => \$rev,
  'maximum_alignment_length=s' => \$maximum_alignment_length,
  'total_mismatches_cutoff=i'  => \$total_mismatches_cutoff,
  'keep_target_suffix'         => \$keep_target_suffix,
  'verbose'                    => \$verbose,
  'help'                       => \$help
);
if ($help) { print $usage; exit; }
$E_cutoff = 5.0 unless($E_cutoff);
$num_mismatch_seed = 2 unless($num_mismatch_seed);
$hsp_cutoff = 17 unless($hsp_cutoff);
$gap_cutoff = 1 unless($gap_cutoff);
$total_mismatches_cutoff = 8 unless($total_mismatches_cutoff);
$GUs_cutoff = 6 unless($GUs_cutoff);
$maximum_alignment_length = 21 unless($maximum_alignment_length);
die("--infile missing...\n") unless($infile);
$penalty_multiplier = 1.5 unless($penalty_multiplier);

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
my $i = 0;
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
        $hash{$curr_query}{length} = $curr_query_length;

        if($verbose){
            print STDERR "curr_query: ".$curr_query."\n";
        }

        next;
    }
    if($_ =~ m/^>>(\S+) .*\((\d+) nt\)$/){
        $i++;
        $curr_target = $1;
        $hash{$curr_query}{$curr_target."_".$i}{name} = $curr_target;
        $hash{$curr_query}{$curr_target."_".$i}{subject_length} = $2;
        next;
    }
    #Smith-Waterman score: 279; 89.5% identity (100.0% similar) in 19 nt overlap (2-20:52164-52182)
    if($_ =~ m/^Smith-Waterman.*\((\d+)-(\d+):(\d+)-(\d+)\)$/){
        my $diff_start = $1 - 1;
        my $diff_end = $hash{$curr_query}{length} - $2;
        my $start = $3 - $diff_start;
        my $end = $4 + $diff_end;
        my $hsp = $4 - $3;
        $hash{$curr_query}{$curr_target."_".$i}{hsp} = $hsp;
        $hash{$curr_query}{$curr_target."_".$i}{start} = $start;
        $hash{$curr_query}{$curr_target."_".$i}{end} = $end;
        if($rev){
            $hash{$curr_query}{$curr_target."_".$i}{q_start} = $1;
            $hash{$curr_query}{$curr_target."_".$i}{q_end} = $2;
            $hash{$curr_query}{$curr_target."_".$i}{strand} = "-";
        }else{
            # If fwd, qstart is in reverse orientation.
            $hash{$curr_query}{$curr_target."_".$i}{q_start} = $2;
            $hash{$curr_query}{$curr_target."_".$i}{q_end} = $1;
            $hash{$curr_query}{$curr_target."_".$i}{strand} = "+";
        }
        next;
    }
    if($_ =~ m/^>--$/){ #same contig, but in another location
        $i++;
        $hash{$curr_query}{$curr_target."_".$i}{name} = $curr_target;
        next;
    }
    if($_ =~ m/^\S+\s+([ACGTU-]*)\s*$/ && $counter_match_nucl_string == 0){
        # Found a query nucl alignment string.
        $curr_start = $-[1];
        $curr_end = $+[1];
        my $str = substr($_, $curr_start, $curr_end - $curr_start);
        $hash{$curr_query}{$curr_target."_".$i}{query_aln} = $str;
        
        $counter_match_nucl_string = 1;
        next;
    }
    if($_ =~ m/^\s+[\.\:]/ && $counter_match_nucl_string == 1){
        # Found a match string.
        #extract substring at previously found positions.
        my $str = substr($_, $curr_start, $curr_end - $curr_start);
        $hash{$curr_query}{$curr_target."_".$i}{match_aln} = $str;
        $str =~ s/^\s+//;
        $str =~ s/\s+$//;
        next;
    }
    if($_ =~ m/^\S+\s+[ACGTBHUYNWRMKS-]*\s*$/ && $counter_match_nucl_string == 1){
        #extract substring at previously found positions.
        my $str = substr($_, $curr_start, $curr_end - $curr_start);
        # Found subject string of current match.
        $hash{$curr_query}{$curr_target."_".$i}{subject_aln} = $str;
       
        if(!exists($hash{$curr_query}{$curr_target."_".$i}{query_aln})){
            die "Problem at $curr_query\n$curr_target"."_".$i."\n".$_."\n"."line number ".$.." in the file\n";
        }

        $counter_match_nucl_string = 0;
        next;
    }
}
close(IN);

print STDERR "Done parsing ssearch standard format outfmt file. Will now parse each alignment...\n";

# Finally sort and print hash.
foreach my $query (keys %hash) {

    foreach my $subject (keys %{ $hash{$query} }) {

        next if($subject eq "length");
        my $aln_str; my $query_str; my $subject_str; my $seed_region_start; my $seed_region_end;
        my $aln_str2; my $query_str2; my $subject_str2;
        
        print STDERR "Processing: ".$query."\t".$subject."\n" if($verbose);

        my $hsp_config = 1; my $offset; my $end;
        if($hsp_config == 1){
            $offset = 0;
            $end = 19;
        }elsif($hsp_config == 2){
            $offset = 1;
            $end = 20;
        }elsif($hsp_config == 3){
            $offset = 2;
            $end = 20;
        }elsif($hsp_config == 4){
            $offset = length($hash{$query}{$subject}{query_aln}) - $hash{$query}{$subject}{q_start};
            $end = 19;
        }

        # If alignments (5'-miRNA-3') vs 5'-DNA-3' are being processed, no need to reverse strings as they already are in the correct orientation.
        if($rev){
            $aln_str = $hash{$query}{$subject}{match_aln};
            $query_str = $hash{$query}{$subject}{query_aln};
            $subject_str = $hash{$query}{$subject}{subject_aln};
            
            $aln_str2       = substr($aln_str,     $offset, $end);
            $query_str2     = substr($query_str,   $offset, $end);
            $subject_str2   = substr($subject_str, $offset, $end);
        
            #$aln_str2 = substr($aln_str, $hash{$query}{$subject}{q_start}, $hash{$query}{$subject}{q_end});
            #$query_str2 = substr($aln_str, $hash{$query}{$subject}{q_start}, $hash{$query}{$subject}{q_end});
            #$subject_str2 = substr($aln_str, $hash{$query}{$subject}{q_start}, $hash{$query}{$subject}{q_end});
        
        # If alignments (3'-miRNA-5') vs 5'-DNA-3' are being processed, we have to reverse strings to make sure that we are starting from the 5' end.
        }else{
            print STDERR "     orig aln string:    ".$hash{$query}{$subject}{match_aln}."\n" if($verbose);
            print STDERR "     query aln string:   ".$hash{$query}{$subject}{query_aln}."\n" if($verbose);
            print STDERR "     subject aln string: ".$hash{$query}{$subject}{subject_aln}."\n" if($verbose);
          
            $aln_str       = $hash{$query}{$subject}{match_aln};
            $query_str     = $hash{$query}{$subject}{query_aln};  
            $subject_str   = $hash{$query}{$subject}{subject_aln};
            $aln_str2      = reverse($aln_str);
            $query_str2    = reverse($query_str);
            $subject_str2  = reverse($subject_str);
            $aln_str2       = substr($aln_str2,     $offset, $end);
            $query_str2     = substr($query_str2,   $offset, $end);
            $subject_str2   = substr($subject_str2, $offset, $end);
        }
        print STDERR "    length aln_str:      ".length($aln_str)."\n" if($verbose);
        print STDERR "    length query_str:    ".length($query_str)."\n" if($verbose);
        print STDERR "    length subject_str:  ".length($subject_str)."\n" if($verbose);
        print STDERR "    aln_str   :          ".$aln_str."\n" if($verbose);
        print STDERR "    query_str :          ".$query_str."\n" if($verbose);
        print STDERR "    aln_str   :          ".$subject_str."\n" if($verbose);
        print STDERR "    length aln_str2:     ".length($aln_str2)."\n" if($verbose);
        print STDERR "    length query_str2:   ".length($query_str2)."\n" if($verbose);
        print STDERR "    length subject_str2: ".length($subject_str2)."\n" if($verbose);
        print STDERR "    aln_str2   :         ".$aln_str2."\n" if($verbose);
        print STDERR "    query_str2 :         ".$query_str2."\n" if($verbose);
        print STDERR "    aln_str2   :         ".$subject_str2."\n" if($verbose);
        print STDERR "------------------------------------------\n" if($verbose); 

        # When looping through hsp, we have to keep the aln string part that matches to the actual hsp
        my @query_char = split(//, $query_str2);
        my @subject_char = split(//, $subject_str2);
        my @aln_char = split(//, $aln_str2);

        my $length = scalar(@aln_char);

        # Compute score according fo Fahlgren 2009 (and psRNATarget).
        # Careful here because in psRNATarget, what they refer to as G:U pairs actually corresponds
        # to . characters in the alignment strings given by SSEARCH36. Furthermore, G:U pairs (or . pairs) 
        # are NOT multiplied by the penalty error multiplier (default 1.5) in the seed region (2-13 nt), only mismatches are multipled.
        # The following routine gives the exact same penalty scoring results
        # to the ones given by psRNATarget (which they call E or Expect value).
        my $z = 1;
        my $score = 0;
        my $gap_status_q = "closed";
        my $gap_status_s = "closed";
        my $gap_status = "closed";
        my $seed_mismatches = 0;
        my $total_mismatches = 0;
        my $gaps = 0;
        my $GUs = 0;
        my $gap_penalty = 0;
        print STDERR "    ".$hash{$query}{$subject}{q_start}."\n" if($verbose);
        for my $el (@aln_char){
            #if($z >= $hash{$query}{$subject}{q_start}){ last; }
            my $el2 = shift(@query_char);
            my $el3 = shift(@subject_char);

            my $is_GU = "no";
            my $curr_score = 0;
            my $mismatch_penalty = 0;
            my $is_mismatch = "no";
            
            if($el eq " " && $el2 ne "-" && $el3 ne "-"){
                $curr_score = 1;
                $is_mismatch = "yes";
                $total_mismatches++;
            }elsif($el eq " "  && $el2 eq "-" && $gap_status_q eq "closed"){ #gap opening
                $curr_score = 2;
                $gap_status_q = "open";
                $is_mismatch = "yes";
                $gaps++;
                $total_mismatches++;
            }elsif($el eq " "  && $el3 eq "-" && $gap_status_s eq "closed"){ #gap opening
                $curr_score = 2;
                $gap_status_s = "open";
                $is_mismatch = "yes";
                $gaps++;
                $total_mismatches++;
            }elsif($el eq " "  && $el2 eq "-" && $gap_status_q eq "open"){ #gap extension
                $curr_score = 0.5;
                $is_mismatch = "yes";
                $gaps++;
                $total_mismatches++;
            }elsif($el eq " "  && $el3 eq "-" && $gap_status_s eq "open"){ #gap extension
                $curr_score = 0.5;
                $is_mismatch = "yes";
                $gaps++;
                $total_mismatches++;
            }elsif($el eq "."){
                $curr_score = 0.5;
                $GUs++;
                $is_GU = "yes";
                $is_mismatch = "no";
                $gap_status_q = "closed";
                $gap_status_s = "closed";
            }elsif($el eq ":"){
                $gap_status_q = "closed";
                $gap_status_s = "closed";
            }else{
                print STDERR "    Did not match anything....: ".$el." ".$el2." ".$el3."\n" if($verbose);
                die;
            }
            
            if($z >= 2 && $z <= 13){
                if($is_mismatch eq "yes"){
                    print STDERR ("    is mismatch in the seed region - will apply penalty... curr_score: ".$curr_score." * 1.5 => ".($curr_score*$penalty_multiplier)."\n") if($verbose);
                    $curr_score = $curr_score * $penalty_multiplier;
                    $seed_mismatches++;
                }else{
                    $curr_score = $curr_score * 1.0;
                }
            }
            $score = $score + $curr_score;
            print STDERR "    $z:$curr_score    $score    $el2 $el $el3\n" if ($verbose);
            $z++;
        }
        $hash{$query}{$subject}{penalty_score} = $score;
        $hash{$query}{$subject}{num_seed_mismatch} = $seed_mismatches;
        $hash{$query}{$subject}{gaps} = $gaps;
        $hash{$query}{$subject}{total_mismatches} = $total_mismatches;
        $hash{$query}{$subject}{GUs} = $GUs;
    }
}

# Print results to standard output
print STDOUT "#query_id\tsubject_id\tmatch_aln\tquery_aln\tsubject_aln\tq_start\tq_end\ts_start\ts_end\texpect_value\tstrand\n";
my %stats;
foreach my $query (keys %hash) {

    foreach my $subject (keys %{ $hash{$query} }) {
        my $subject_id = $subject;
        unless($keep_target_suffix){
            $subject_id =~ s/_\d+$//;
        }

        next if($subject eq "length");
        if($hash{$query}{$subject}{penalty_score} <= $E_cutoff && 
            $hash{$query}{$subject}{num_seed_mismatch} <= $num_mismatch_seed && 
            $hash{$query}{$subject}{hsp} >= $hsp_cutoff &&
            $hash{$query}{$subject}{gaps} <= $gap_cutoff &&
            $hash{$query}{$subject}{total_mismatches} <= $total_mismatches_cutoff &&
            length($hash{$query}{$subject}{subject_aln}) <= $maximum_alignment_length){
            
            if($hash{$query}{$subject}{GUs} <= $GUs_cutoff){
                print STDOUT $query."\t".$subject_id."\t".$hash{$query}{$subject}{match_aln}."\t".$hash{$query}{$subject}{query_aln}."\t".$hash{$query}{$subject}{subject_aln}."\t".$hash{$query}{$subject}{q_start}."\t".$hash{$query}{$subject}{q_end}."\t".$hash{$query}{$subject}{start}."\t".$hash{$query}{$subject}{end}."\t".$hash{$query}{$subject}{penalty_score}."\t".$hash{$query}{$subject}{strand}."\n";
            }else{
                $stats{GU_reject}++;
            }
        }else{
            $stats{other_reject}++;
               print STDERR $query."\t".$subject_id."\t".$hash{$query}{$subject}{match_aln}."\t".$hash{$query}{$subject}{query_aln}."\t".$hash{$query}{$subject}{subject_aln}."\t".$hash{$query}{$subject}{q_start}."\t".$hash{$query}{$subject}{q_end}."\t".$hash{$query}{$subject}{start}."\t".$hash{$query}{$subject}{end}."\t".$hash{$query}{$subject}{penalty_score}."\t".$hash{$query}{$subject}{strand}."\n";
        }
    }
}
print STDERR Dumper(\%stats) if($verbose);
if($verbose){ print STDERR  Dumper(\%hash); }

sub complement_IUPAC {
   my $dna = shift;

   # complement the reversed DNA sequence
   my $comp = $dna;
   $comp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy-/TVGHCDKNYSAABWXRtvghcdknysaabwxr-/;
   return $comp;
}

exit;


