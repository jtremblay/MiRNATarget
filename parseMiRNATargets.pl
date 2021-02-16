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
From the default output of ssearch36, distribute miRNA into their site types. According to 
scheme found in PMC3499661. Results of this script matches almost exactly the ones outputted 
by psRNATarget. The scoring sheme gives exactly the same results as psRNATarget.

INPUT:
--infile <string>               : ssearch default out format file
--E_cutoff <float>              : default = 5.0  - Expectation value is the penalty for the mismatches 
                                  between mature small RNA and the target sequence. A higher value indicates 
                                  less similarity (and possibility) between small RNA and the target candidate. 
                                  The default penalty rule is set up by the scoring schema. Maximum expectation is the cutoff; 
                                  any small RNA-target pair with an expectation less than the cutoff will be discarded in the final result. 
                                  The recommended values are 3.0-5.0 depending on the scoring schema. 
--num_mismatch_seed <int>       : default = 2 - Maximum of allowed mismatches excluding G:U pairs.
--hsp_cutoff <int>              : default = 17 - HSPs shorter than this value will be discarded.
--gap_cutoff <int>              : default = 1 - alignments having more than <int> gaps willbe discarded.
--total_mismatches_cutoff <int> : default = 8 - alignments showing more than <int> mismatches will be discarded.
--GUs_cutoff <int>              : default = 6 - alignments having more than <int> mismatches will be discarded.
--keep_target_suffix            : set flag if you wish to keep the temporary _<int> suffix appended at the end
                                  of subjects IDs.
--rev                           : Set flag if alignments of miRNAs were done on revcomp subject sequences.

OUTPUT:
STDOUT <string>    : standard output

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@nrc-cnrc.gc.ca

ENDHERE

## OPTIONS
my ($help, $infile, $E_cutoff, $num_mismatch_seed, $hsp_cutoff, $gap_cutoff, $total_mismatches_cutoff, $GUs_cutoff, $keep_target_suffix, $rev);
my $verbose = 0;

GetOptions(
  'infile=s'                  => \$infile,
  'E_cutoff=f'                => \$E_cutoff,
  'num_mismatch_seed=i'       => \$num_mismatch_seed,
  'hsp_cutoff=i'              => \$hsp_cutoff,
  'gap_cutoff=i'              => \$gap_cutoff,
  'GUs_cutoff=i'              => \$GUs_cutoff,
  'rev'                       => \$rev,
  'total_mismatches_cutoff=i' => \$total_mismatches_cutoff,
  'keep_target_suffix'        => \$keep_target_suffix,
  'verbose'                   => \$verbose,
  'help'                      => \$help
);
if ($help) { print $usage; exit; }
$E_cutoff = 5.0 unless($E_cutoff);
$num_mismatch_seed = 2 unless($num_mismatch_seed);
$hsp_cutoff = 17 unless($hsp_cutoff);
$gap_cutoff = 1 unless($gap_cutoff);
$total_mismatches_cutoff = 8 unless($total_mismatches_cutoff);
$GUs_cutoff = 6 unless($GUs_cutoff);
die("--infile missing...\n") unless($infile);

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
        $counter_match_nucl_string = 0;
        $curr_start = 0;
        $curr_end = 0;
        $curr_query = $1;
        $curr_query_length = $2;
        $i = 0;
        $hash{$curr_query}{length} = $curr_query_length;
        next;
    }
    if($_ =~ m/^>>(\S+) .*\((\d+) nt\)$/){
        $i++;
        $curr_target = $1;
        $hash{$curr_query}{$curr_target."_".$i}{name} = $curr_target;
        $hash{$curr_query}{$curr_target."_".$i}{subject_length} = $2;
        next;
    }
    if($_ =~ m/^>--$/){
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
		#my @str = split(//, $str);
		$str =~ s/^\s+//;
		$str =~ s/\s+$//;
		my $hsp = length($str);
		$hash{$curr_query}{$curr_target."_".$i}{hsp} = $hsp;
        next;
    }
    #Smith-Waterman score: 279; 89.5% identity (100.0% similar) in 19 nt overlap (2-20:52164-52182)
    if($_ =~ m/^Smith-Waterman.*\((\d+)-(\d+):(\d+)-(\d+)\)$/){
		my $diff_start = $1 - 1;
		my $diff_end = $hash{$curr_query}{length} - $2;
        my $start = $3 - $diff_start;
        my $end = $4 + $diff_end;
		my $hsp = $4 - $3;
        #print STDERR "1: ".$start."\n";
        #print STDERR "2: ".$end."\n";
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
    if($_ =~ m/^\S+\s+[ACGTUYNWRMKS-]*\s*$/ && $counter_match_nucl_string == 1){
        #extract substring at previously found positions.
        my $str = substr($_, $curr_start, $curr_end - $curr_start);
        # Found subject string of current match.
        $hash{$curr_query}{$curr_target."_".$i}{subject_aln} = $str;
       
        if(!exists($hash{$curr_query}{$curr_target."_".$i}{query_aln})){
            print STDERR "line # ".$.."\n" if($verbose);
            print STDERR Dumper(\%hash) if($verbose);
            die "Problem at $curr_query\n$curr_target"."_".$i."\n";
        }

        $counter_match_nucl_string = 0;
        next;
    }
}
close(IN);

print STDERR "Done parsing ssearch standard format outfmt file. Will now parse each alignments...\n";

# Finally sort and print hash.
#my %hash_aln;
foreach my $query (keys %hash) {

    foreach my $subject (keys %{ $hash{$query} }) {

        next if($subject eq "length");
        my $aln_str; my $query_str; my $subject_str;

        # If alignments (5'-miRNA-3') vs 5'-DNA-3' are being processed, no need to reverse strings as they already are in the correct orientation.
        if($rev){
            $aln_str = $hash{$query}{$subject}{match_aln};
            $query_str = $hash{$query}{$subject}{query_aln};
            $subject_str = $hash{$query}{$subject}{subject_aln};
        }else{
            $aln_str = reverse($hash{$query}{$subject}{match_aln});
            $query_str = reverse($hash{$query}{$subject}{query_aln});
            $subject_str = reverse($hash{$query}{$subject}{subject_aln});
        }
        my @query_char = split(//, $query_str);
        my @subject_char = split(//, $subject_str);
		my @aln_char = split(//, $aln_str);

        my $length = scalar(@aln_char);

        #$hash_aln{reverse($aln_str)}++;
        
		# Compute score according fo Fahlgren 2009 (and psRNATarget).
 		# The following routine gives the exact same penalty scoring results
        # to the ones given by psRNATarget (which they call E or Expect value.
		print STDERR $query."\t".$subject."\n" if($verbose);
		
		my $z = 1;
		my $score = 0;
        my $gap_status_q = "closed";
        my $gap_status_s = "closed";
        my $seed_mismatches = 0;
        my $total_mismatches = 0;
		my $gaps = 0;
        my $GUs = 0;
		
        for my $el (@aln_char){
            my $el2 = shift(@query_char);
            my $el3 = shift(@subject_char);

			my $is_GU = "no";
			my $curr_score = 0;
            my $is_mismatch = "no";
			
			# Implement $el3
			if($el eq " " && $el2 ne "-" && $el3 ne "-" && $z == $length){
				$curr_score = 0;
            }elsif($el eq " " && $el2 ne "-" && $el3 ne "-"){
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
                $curr_score = 1.5;
                $is_mismatch = "yes";
                $gaps++;
                $total_mismatches++;
            }elsif($el eq " "  && $el3 eq "-" && $gap_status_s eq "open"){ #gap extension
                $curr_score = 1.5;
                $is_mismatch = "yes";
                $gaps++;
                $total_mismatches++;
            }elsif($el eq "."){
				$curr_score = 0.5;
                $is_mismatch = "no";
                $total_mismatches++;
				$GUs++;
                $gap_status_q = "closed";
                $gap_status_s = "closed";
				$is_GU = "yes";
			}elsif($el eq ":"){
                $gap_status_q = "closed";
                $gap_status_s = "closed";
            }
			
            # adjust if we are in seed region.
            if($z >= 2 && $z <= 13){
				if($is_GU eq "yes"){
					$curr_score = $curr_score * 1.0;
				}elsif($is_GU eq "no"){
					$curr_score = $curr_score * 1.5;
				}
                if($is_mismatch eq "yes"){
                    $seed_mismatches++;
                }
			}
			$score = $score + $curr_score;
		    print STDERR "$z:$curr_score    $score\n" if ($verbose);
			$z++;
        }
		$hash{$query}{$subject}{penalty_score} = $score;
		$hash{$query}{$subject}{num_seed_mismatch} = $seed_mismatches;
		$hash{$query}{$subject}{gaps} = $gaps;
		$hash{$query}{$subject}{total_mismatches} = $total_mismatches;
		$hash{$query}{$subject}{GUs} = $GUs;
    }
}
print STDERR Dumper(\%hash) if($verbose);

# Take 20 most abundant matches.
#my %hash_aln_ma;
#foreach my $match_string (sort { $hash_aln{$b} <=> $hash_aln{$a} } keys %hash_aln) {
#    #next if($i > 200);
#    $hash_aln_ma{$match_string}{freq} = $hash_aln{$match_string};
#    print STDERR $match_string."=>".$hash_aln{$match_string}."\n" if($verbose);
#    #$i++;
#}
#print STDERR Dumper(\%hash_aln_ma);
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
            $hash{$query}{$subject}{total_mismatches} <= $total_mismatches_cutoff){
			
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

sub complement_IUPAC {
   my $dna = shift;

   # complement the reversed DNA sequence
   my $comp = $dna;
   $comp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy-/TVGHCDKNYSAABWXRtvghcdknysaabwxr-/;
   return $comp;
}

exit;


### If ssearch output 9C:
#while(<IN>){
#    chomp;
#    if($_ =~ m/^#/){
#        next;
#    }
#    if($_ =~ m/^>>><<</){
#        $counter_match_nucl_string = 0;
#        $curr_start = 0;
#        $curr_end = 0;
#        next;
#    }
#    if($_ =~ m/^>>>(\S+), (\d+) nt/){
#        $curr_query = $1;
#        $curr_query_length = $2;
#        $hash{$curr_query}{length} = $curr_query_length;
#        next;
#    }
#    if($_ =~ m/^>>(\S+) .*\(\d+ nt\)$/){
#        $curr_target = $1;
#        #$curr_target =~ s/\s+$//;
#        #print STDERR $curr_target."\n";
#        $hash{$curr_query}{$curr_target}{name} = $curr_target;
#        next;
#    }
#    #GL635794.1 Neisseria mucosa C102 genomic scaffold superc (494621) [f]
#    if($_ =~ m/^\S+\s.*\(\d+\) \[/){
#        #print STDERR $_."\n";
#        next;
#    }
#    if($_ =~ m/^\S+\s+([ACGTU-]*)\s*$/ && $counter_match_nucl_string == 0){
#        # Found a query nucl alignment string.
#        $curr_start = $-[1];
#        $curr_end = $+[1];
#        #print STDERR "curr_start: ".$curr_start."\n";
#        #print STDERR "curr_end: ".$curr_end."\n";
#        my $str = substr($_, $curr_start, $curr_end - $curr_start);
#        $hash{$curr_query}{$curr_target}{query_aln} = $str;
#        
#        $counter_match_nucl_string = 1;
#        next;
#    }
#    if($_ =~ m/^\s+[\.\:]/ && $counter_match_nucl_string == 1){
#        # Found a match string.
#        #extract substring at previously found positions.
#        my $str = substr($_, $curr_start, $curr_end - $curr_start);
#        $hash{$curr_query}{$curr_target}{match_aln} = $str;
#        next;
#    }
#    #Smith-Waterman score: 279; 89.5% identity (100.0% similar) in 19 nt overlap (2-20:52164-52182)
#    if($_ =~ m/^Smith-Waterman.*\(\d+-\d+:(\d+)-(\d+)\)$/){
#        my $start = $1;
#        my $end = $2;
#        #print STDERR "1: ".$start."\n";
#        #print STDERR "2: ".$end."\n";
#        $hash{$curr_query}{$curr_target}{start} = $start;
#        $hash{$curr_query}{$curr_target}{end} = $end;
#        next;
#    }
#    if($_ =~ m/^\S+\s+[ACGTUYNWRMKS-]*\s*$/ && $counter_match_nucl_string == 1){
#        #extract substring at previously found positions.
#        my $str = substr($_, $curr_start, $curr_end - $curr_start);
#        # Found subject string of current match.
#        $hash{$curr_query}{$curr_target}{subject_aln} = $str;
#       
#        if(!exists($hash{$curr_query}{$curr_target}{query_aln})){
#            print STDERR "line # ".$.."\n";
#            print STDERR Dumper(\%hash);
#            die "Problem at $curr_query\n$curr_target\n";
#        }
#
#        $counter_match_nucl_string = 0;
#        next;
#    }
#}
#close(IN);
