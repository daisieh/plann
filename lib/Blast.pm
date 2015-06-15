#!/usr/bin/perl

package Blast;
use strict;
use FindBin;
use lib "$FindBin::Bin/..";
use Subfunctions qw(split_seq reverse_complement meld_matrices);
use Data::Dumper;
use XML::Simple;
#no, we should use XML::XPath!


BEGIN {
	require Exporter;
	# set the version for version checking
	our $VERSION     = 1.00;
	# Inherit from Exporter to export functions and variables
	our @ISA         = qw(Exporter);
	# Functions and variables which are exported by default
	our @EXPORT      = qw(blast_to_ref debug);
	# Functions and variables which can be optionally exported
	our @EXPORT_OK   = qw(parse_xml revcomp_hsp sort_hsps_by_score sort_regions_by_start sort_hsps_by_hit_start);
}
my $debug = 0;

sub debug {
	my $str = shift;

	if ($debug == 1) {
		print STDERR "$str";
	}
}


sub blast_to_ref {
	my $blast_xml = shift;
	my %result_matrix = ();

	my $hit_array = parse_xml($blast_xml);

	foreach my $hit (@$hit_array) {
		my $query_name = $hit->{'query'}->{'name'};
		my $hsps = $hit->{'hsps'}; # Key Hsp has a value that is an anonymous array of the hit hashes.
		my $subject_name = $hit->{'subject'}->{'name'};
		$subject_name =~ s/\s+/_/g;
		debug ("HIT $subject_name QUERY $query_name\n");
		my @sorted_hsps = sort { $b->{'score'} - $a->{'score'} } @$hsps;
		my @selected_hsps = ();
		my ($query_start, $query_end, $hit_start, $hit_end, $align_strand) = 0;
		while (my $hsp = shift @sorted_hsps) {
			# if $hsp is on the minus strand, reverse-comp before dealing with it.
			if ($hsp->{'hit-frame'} < 0) {
				revcomp_hsp($hsp);
			}
			my $aln_length = $hsp->{'align-len'};
			my $aln_percent = sprintf("%.2f",$hsp->{'identity'} / $aln_length);
	# 					if ($hsp->{'evalue'} > 10) {
	# 						# this is no good: move on to the next hit.
	# 						next;
	# 					}
			if (@selected_hsps == 0) {
				# this is the first chunk of sequence.
				$query_start = $hsp->{'query-from'};
				$query_end = $hsp->{'query-to'};
				$hit_start = $hsp->{'hit-from'};
				$hit_end = $hsp->{'hit-to'};
				$align_strand = $hsp->{'hit-frame'};
				debug ("\tfirst seq with score ".$hsp->{'score'}.", strand $align_strand: ".$hsp->{'hit-from'} ."-".$hsp->{'hit-to'}."\t".$hsp->{'query-from'}."-".$hsp->{'query-to'}."\n");
				push @selected_hsps, $hsp;
			} else {
				# if it's not the first chunk, we need to check to see if it falls in a reasonable part of the sequence: both hit and query have to lie before or after the working region.
				debug ("looking for hsps that work with $align_strand: $hit_start-$hit_end\t$query_start-$query_end\n");
				debug ("\tlooking at next seq with score ".$hsp->{'score'}.", strand ".$hsp->{'hit-frame'}.": ".$hsp->{'hit-from'} ."-".$hsp->{'hit-to'}."\t".$hsp->{'query-from'}."-".$hsp->{'query-to'}."\n");
				if ($hsp->{'hit-frame'} != $align_strand) {
					debug ("\t\tNO: next seq on wrong strand\n");
				} else {
					if ($align_strand > 0) {
						if (($hsp->{'hit-to'} < $hit_start) && ($hsp->{'query-to'} < $query_start)) {
							# we know $hsp will fall completely on the left, so we're good.
							debug ("\t\tYES: next seq falls to the left of $hit_start and $query_start: ".$hsp->{'hit-from'} ."-".$hsp->{'hit-to'}."\t".$hsp->{'query-from'}."-".$hsp->{'query-to'}."\n");
							$query_start = $hsp->{'query-from'};
							$hit_start = $hsp->{'hit-from'};
							push @selected_hsps, $hsp;
						} elsif (($hsp->{'hit-from'} > $hit_end) && ($hsp->{'query-from'} > $query_end)) {
							# we know $hsp will fall completely on the left, so we're good.
							debug ("\t\tYES: next seq falls to the right of $hit_end and $query_end: ".$hsp->{'hit-from'} ."-".$hsp->{'hit-to'}."\t".$hsp->{'query-from'}."-".$hsp->{'query-to'}."\n");
							$query_end = $hsp->{'query-to'};
							$hit_end = $hsp->{'hit-to'};
							push @selected_hsps, $hsp;
						} else {
							debug ("\t\tNO: hit and query fell on opposite sides\n");
						}
					} else {
						# we're dealing with a minus strand:
						# hit should be the same, but query is backwards.
						if (($hsp->{'hit-to'} < $hit_start) && ($hsp->{'query-from'} > $query_end)) {
							# we know $hsp will fall completely on the left, so we're good.
							debug ("\t\tYES: next seq falls to the left of $hit_start and -$query_end: ".$hsp->{'hit-from'} ."-".$hsp->{'hit-to'}."\t".$hsp->{'query-from'}."-".$hsp->{'query-to'}."\n");
							$query_end = $hsp->{'query-to'};
							$hit_start = $hsp->{'hit-from'};
							push @selected_hsps, $hsp;
						} elsif (($hsp->{'hit-from'} > $hit_end) && ($hsp->{'query-to'} < $query_start)) {
							# we know $hsp will fall completely on the left, so we're good.
							debug ("\t\tYES: next seq falls to the right of $hit_end and -$query_start: ".$hsp->{'hit-from'} ."-".$hsp->{'hit-to'}."\t".$hsp->{'query-from'}."-".$hsp->{'query-to'}."\n");
							$query_start = $hsp->{'query-from'};
							$hit_end = $hsp->{'hit-to'};
							push @selected_hsps, $hsp;
						} else {
							debug ("\t\tNO: hit and query fell on opposite sides\n");
						}
					}
				}
			}

			# remove all hits that would fall within this range, starting from the far end of the array (this ensures we don't mess up the count when we splice things out).
			my $this_hsp = 0;
			debug ("\t\tremoving any hsps that fall in the current range: $hit_start-$hit_end\t$query_start-$query_end\n");
			for (my $j=@sorted_hsps-1; $j>=0; $j--) {
				$this_hsp = $sorted_hsps[$j];
				if (($this_hsp->{'hit-from'} >= $hit_start) && ($this_hsp->{'hit-to'} <= $hit_end)) {
					debug ("\t\t\tremoving hsp ".$this_hsp->{'hit-from'}."-".$this_hsp->{'hit-to'}." with score ". $this_hsp->{'score'} . " because it's inside the hit\n");
					splice (@sorted_hsps,$j,1);
				} elsif (($this_hsp->{'query-from'} >= $query_start) && ($this_hsp->{'query-to'} <= $query_end)) {
					debug ("\t\t\tremoving hsp ".$this_hsp->{'query-from'}."-".$this_hsp->{'query-to'}." with score ". $this_hsp->{'score'} . " because it's inside the query\n");
					splice (@sorted_hsps,$j,1);
				} else {
					debug ("\t\t\thsp ".$this_hsp->{'hit-from'}."-".$this_hsp->{'hit-to'}."\t".$this_hsp->{'query-from'}."-".$this_hsp->{'query-to'}." is not inside\n");
				}
			}
			debug ("\tthere are still ". @sorted_hsps . " hsps left\n");
		}
		debug ("final sequence:\n");
		my $hit_end = 0;
		my $sequence = "";
		my @selected_hsps = sort { $a->{'hit-from'} - $b->{'hit-from'} } @selected_hsps;
		my $total_identity = 0;
		my $total_aln_length = 0;
		foreach my $hsp (@selected_hsps) { # for each hsp for this query
			my $aln_length = $hsp->{'align-len'};
			my $aln_percent = sprintf("%.2f",$hsp->{'identity'} / $aln_length);
			debug ("\t".$hsp->{'hit-from'}."-".$hsp->{'hit-to'}."\t".$hsp->{'query-from'}."-".$hsp->{'query-to'}."\t$aln_length\t$aln_percent\t".$hsp->{'score'}."\n");
			$total_aln_length += $aln_length;
			$total_identity += $hsp->{'identity'};
			my $last_end = $hit_end;
			$hit_end = $hsp->{'hit-to'};
			# remove gaps from query seq that might have been inserted into the ref seq.
			while ($hsp->{'hseq'} =~ /^(.*?)(-+)(.*)$/) {
				my $left = length($1);
				my $gap = length($2);
				my $right = length($3);
				$hsp->{'hseq'} = $1 . $3;
# 				$hsp->{'qseq'} =~ /^(.{$left})(.{$gap})(.{$right})/;
				my ($startseq, $regionseq, $endseq) = split_seq ( $hsp->{'qseq'}, $left+1, $left+$gap);
				$hsp->{'qseq'} = $startseq . $endseq;
			}
			$sequence .= "n" x ($hsp->{'hit-from'} - $last_end - 1) . $hsp->{'qseq'};
		}
		debug ("aligned $total_identity out of $total_aln_length chars\n");
		if ($sequence =~ /\w/) {
			$result_matrix{$subject_name}->{$query_name} = $sequence;
			debug ("put seq in result_matrix at $subject_name->$query_name\n");
		} else {
			debug ("\tno match\n");
		}
	}
	return \%result_matrix;
}

# parse_xml takes the freakshow that is the BLAST XML format and simplifies it:
# the result is an array of hits:
# each hit is a hash with three keys: query, subject, hsps.
# the "hsps" value is an array of high-scoring segment pairs (hsps).
# each hsp has a bunch of associated key/value pairs.

sub parse_xml {
	my $blast_xml = shift;

	my $parser = new XML::Simple();

	open XML_FH, "<:crlf", $blast_xml;

	my $xml = "";
	my $query_name = "";

	# let's just pretend that there could never be more than one xml string.
	while (my $line = readline XML_FH) {
		$xml .= $line;
	}

	close XML_FH;
	my @hit_array = ();
	my $hit_hash = {};
	my $tree = $parser->XMLin($xml, ForceArray => 1);
	my $iterations = ($tree->{'BlastOutput_iterations'}[0]->{'Iteration'}); # key Iteration represents a single blast search; has a value that is an anonymous array of iteration hashes.
# 	print @{$tree->{'BlastOutput_iterations'}[0]->{'Iteration'}} . " iterations\n";
	foreach my $iteration (@$iterations) {
		foreach my $hit (@{$iteration->{'Iteration_hits'}[0]->{'Hit'}}) { # each Hit in the iteration can have multiple Hsps.
			my $hit_hash = {};
			push @hit_array, $hit_hash;
			$hit_hash->{'query'}->{'name'} = $iteration->{'Iteration_query-def'}[0];
			$hit_hash->{'query'}->{'length'} = $iteration->{'Iteration_query-len'}[0];
			my @hsp_array = ();
			$hit_hash->{'hsps'} = \@hsp_array;
			$hit_hash->{'subject'}->{'name'} = $hit->{'Hit_def'}[0];
			$hit_hash->{'subject'}->{'length'} = $hit->{'Hit_len'}[0];
			foreach my $hsp (@{$hit->{'Hit_hsps'}[0]->{'Hsp'}}) { # Key Hsp has a value that is an anonymous array of the hsp hashes.
				my $hsp_hash = {};
				push @hsp_array, $hsp_hash;
				$hsp_hash->{'hseq'} = $hsp->{'Hsp_hseq'}[0];
				$hsp_hash->{'hit-frame'} = $hsp->{'Hsp_hit-frame'}[0];
				$hsp_hash->{'bit-score'} = $hsp->{'Hsp_bit-score'}[0];
				$hsp_hash->{'evalue'} = $hsp->{'Hsp_evalue'}[0];
				$hsp_hash->{'qseq'} = $hsp->{'Hsp_qseq'}[0];
				$hsp_hash->{'hit-to'} = $hsp->{'Hsp_hit-to'}[0];
				$hsp_hash->{'identity'} = $hsp->{'Hsp_identity'}[0];
				$hsp_hash->{'hit-from'} = $hsp->{'Hsp_hit-from'}[0];
				$hsp_hash->{'query-from'} = $hsp->{'Hsp_query-from'}[0];
				$hsp_hash->{'gaps'} = $hsp->{'Hsp_gaps'}[0];
				$hsp_hash->{'num'} = $hsp->{'Hsp_num'}[0];
				$hsp_hash->{'query-to'} = $hsp->{'Hsp_query-to'}[0];
				$hsp_hash->{'positive'} = $hsp->{'Hsp_positive'}[0];
				$hsp_hash->{'query-frame'} = $hsp->{'Hsp_query-frame'}[0];
				$hsp_hash->{'score'} = $hsp->{'Hsp_score'}[0];
				$hsp_hash->{'align-len'} = $hsp->{'Hsp_align-len'}[0];
				$hsp_hash->{'midline'} = $hsp->{'Hsp_midline'}[0];

			}
		}
	}
	return \@hit_array;
}

sub revcomp_hsp {
	my $hsp = shift;

	my $hit_to = $hsp->{'hit-to'};
	$hsp->{'hit-to'} = $hsp->{'hit-from'};
	$hsp->{'hit-from'} = $hit_to;
# 	my $query_to = $hsp->{'query-to'};
# 	$hsp->{'query-to'} = $hsp->{'query-from'};
# 	$hsp->{'query-from'} = $query_to;
	$hsp->{'qseq'} = lc(reverse_complement($hsp->{'qseq'}));
	$hsp->{'hseq'} = lc(reverse_complement($hsp->{'hseq'}));
	$hsp->{'hit-frame'} = -1 * $hsp->{'hit-frame'};
	return $hsp;
}

###### SORT FUNCTIONS
sub sort_hsps_by_score ($$) {
	my $a = $_[0];
	my $b = $_[1];
	my $score = $b->{'bit-score'} - $a->{'bit-score'};
	if ($score == 0) {
		my $b_direction = ($b->{'query-to'} - $b->{'query-from'})/($b->{'hit-to'} - $b->{'hit-from'});
		my $a_direction = ($a->{'query-to'} - $a->{'query-from'})/($a->{'hit-to'} - $a->{'hit-from'});
		if ($b_direction > $a_direction) {
			$score = 1;
		} elsif ($a_direction < $b_direction) {
			$score = -1;
		} else {
			$score = 0;
		}
	}
	return $score;
}

sub sort_hsps_by_hit_start ($$) {
	my $a = $_[0];
	my $b = $_[1];
	my $b_start = $b->{'hit-from'};
	my $a_start = $a->{'hit-from'};
	return ($a_start <=> $b_start);
}

sub sort_hsps_by_query_start ($$) {
	my $a = $_[0];
	my $b = $_[1];
	my $b_start = $b->{'query-from'};
	my $a_start = $a->{'query-from'};
	return ($a_start <=> $b_start);
}

sub sort_regions_by_start ($$) {
	my $a = $_[0];
	my $b = $_[1];
	$b =~ /.*\t(\d+)\t(\d+)/;
	my $b_start = $2;
	$a =~ /.*\t(\d+)\t(\d+)/;
	my $a_start = $2;
	return ($a_start <=> $b_start);
}



return 1;

