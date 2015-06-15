#!/usr/bin/perl

# WARNING: this only handles a single genbank file per program.
# TODO: object-orient this.

package Genbank;
use strict;
use FindBin;
use lib "$FindBin::Bin/..";
use Subfunctions qw(split_seq reverse_complement);
use Data::Dumper;


BEGIN {
	require Exporter;
	# set the version for version checking
	our $VERSION     = 1.00;
	# Inherit from Exporter to export functions and variables
	our @ISA         = qw(Exporter);
	# Functions and variables which are exported by default
	our @EXPORT      = qw(parse_genbank sequence_for_interval sequin_feature stringify_feature parse_feature_desc parse_interval parse_qualifiers within_interval write_sequin_tbl feature_table_from_genbank set_sequence get_sequence get_name write_features_as_fasta clone_features simplify_genbank_array);
	# Functions and variables which can be optionally exported
	our @EXPORT_OK   = qw();
}

my $line = "";
my $in_features = 0;
my $in_sequence = 0;
my $feat_desc_string = "";
my $sequence = "";
my $curr_gene = {};
my $gb_name = "";

sub set_sequence {
	$sequence = shift;
}

sub get_sequence {
	return $sequence;
}

sub get_name {
	return $gb_name;
}


sub parse_genbank {
	my $gbfile = shift;

	open FH, "<:crlf", $gbfile;

	my @gene_array = ();
	$line = "";
	$in_features = 0;
	$in_sequence = 0;
	$feat_desc_string = "";
	$sequence = "";
	$curr_gene = {};
	$gb_name = "";

	$line = readline FH;
	while (defined $line) {
		if ($line =~ /^\s+$/) {
			# if the line is blank, skip to end
		} elsif ($line =~ /^\/\//) { # we've finished the file
			last;
		} elsif ($line =~ /^\S/) { # we're looking at a new section
			if ($line =~ /DEFINITION\s+(.*)$/) { # if line has the definition, keep this.
				$gb_name = $1;
			} elsif ($line =~ /FEATURES/) { # if we hit a line that says FEATURES, we're starting the features.
				$in_features = 1;
			} elsif ($line =~ /ORIGIN/) { # if we hit a line that says ORIGIN, we have sequence
				# this is the end of the features: process the last feature and finish.
				parse_feature_desc ($feat_desc_string, \@gene_array);
				$in_features = 0;
				$in_sequence = 1;
			}
		} elsif ($in_sequence == 1) {
			$line =~ s/\d//g;
			$line =~ s/\s//g;
			$sequence .= $line;
			# process lines in the sequence.
		} elsif ($in_features == 1) {
			# three types of lines in features:
			if ($line =~ /^\s{5}(\S.+)$/) {
				# start of a feature descriptor. Parse whatever feature and qualifiers we might've been working on and move on.
				if ($feat_desc_string ne "") {
					my $feat_desc_hash = parse_feature_desc ($feat_desc_string, \@gene_array);
				}
				$feat_desc_string = $1;
			} elsif ($line =~ /^\s{21}(\S.+)$/) {
				$feat_desc_string .= "+$1";
			}
		}

		$line = readline FH;
	}
	close FH;
	return \@gene_array;
}

sub simplify_genbank_array {
	my $gene_array = shift;

	my @simplified_array = ();
	foreach my $feat (@$gene_array) {
		if ($feat->{'type'} =~ /gene|source/) {
			my @newcontains = ();
			foreach my $subfeat (@{$feat->{'contains'}}) {
				if ($subfeat->{'type'} !~ /exon|intron|STS/) {
					push @newcontains, $subfeat;
				}
			}
			$feat->{'contains'} = \@newcontains;
			push @simplified_array, $feat;
		}
	}
	return \@simplified_array;
}

sub parse_feature_sequences {
	my $gene_array = shift;

	my $gene_id = 0;
	foreach my $gene (@$gene_array) {
		if ($gene->{'type'} eq "gene") {
			my $interval_str = flatten_interval ($gene->{'region'});
			my $geneseq = sequence_for_interval ($interval_str);
			my $genename = $gene->{'qualifiers'}->{'gene'};
			foreach my $feat (@{$gene->{'contains'}}) {
				my $feat_id = 0;
				foreach my $reg (@{$feat->{'region'}}) {
					my $strand = "+";
					my ($start, $end) = split (/\.\./, $reg);
					if ($end < $start) {
						$strand = "-";
						my $oldend = $start;
						$start = $end;
						$end = $oldend;
					}
					my $regseq = sequence_for_interval ($reg);
					my $featname = $feat->{'type'};
					$feat_id++;
				}
			}
			my $interval_str = flatten_interval ($gene->{'region'});
			my @gene_interval = ($interval_str);
			$gene_id++;
		}
	}
}

sub clone_features {
	my $gene_array = shift;
	my $flattened_hash = {};
	my @flattened_names = ();
	my $gene_id = 0;

	foreach my $gene (@$gene_array) {
		if ($gene->{'type'} eq "gene") {
			my $interval_str = flatten_interval ($gene->{'region'});
			my @gene_interval = ($interval_str);
			# we need to disambiguate different copies of the gene
			my $genename = $gene->{'qualifiers'}->{'gene'} . "_$gene_id";
			my $strand = "+";
			my $reg = @{$gene->{'region'}}[0];
			my ($start, $end) = split (/\.\./, $reg);
			if ($end < $start) {
				$strand = "-";
				my $oldend = $start;
				$start = $end;
				$end = $oldend;
			}
			my $regseq = sequence_for_interval ($reg);
			push @flattened_names, $genename;
			$flattened_hash->{$genename}->{'strand'} = $strand;
			$flattened_hash->{$genename}->{'start'} = $start;
			$flattened_hash->{$genename}->{'end'} = $end;
			$flattened_hash->{$genename}->{'characters'} = $regseq;
			$flattened_hash->{$genename}->{'gene'} = $genename;
			$flattened_hash->{$genename}->{'type'} = "gene";
			$flattened_hash->{$genename}->{'id'} = $gene_id;
			$flattened_hash->{$genename}->{'qualifiers'} = {};
			add_qualifiers_to_feature_hash($flattened_hash->{$genename}->{'qualifiers'}, $gene);
			$flattened_hash->{$genename}->{'contains'} = [];
			foreach my $feat (@{$gene->{'contains'}}) {
				add_qualifiers_to_feature_hash($flattened_hash->{$genename}->{'qualifiers'}, $feat);
				my $feat_id = 0;
				foreach my $reg (@{$feat->{'region'}}) {
					my $strand = "+";
					my ($start, $end) = split (/\.\./, $reg);
					if ($end < $start) {
						$strand = "-";
						my $oldend = $start;
						$start = $end;
						$end = $oldend;
					}
					my $regseq = sequence_for_interval ($reg);
					my $featname = $feat->{'type'};
					my $fullname = "$genename"."_$featname"."_$feat_id";
					push @flattened_names, $fullname;
					push @{$flattened_hash->{$genename}->{'contains'}}, $fullname;
					$flattened_hash->{$fullname}->{'strand'} = $strand;
					$flattened_hash->{$fullname}->{'start'} = $start;
					$flattened_hash->{$fullname}->{'end'} = $end;
					$flattened_hash->{$fullname}->{'characters'} = $regseq;
					$flattened_hash->{$fullname}->{'gene'} = $genename;
					$flattened_hash->{$fullname}->{'type'} = $featname;
					$flattened_hash->{$fullname}->{'id'} = $gene_id;
					$feat_id++;
				}
			}
			$gene_id++;
		}
	}
	return ($flattened_hash, \@flattened_names);
}

sub write_features_as_fasta {
	my $gene_array = shift;

	my ($flattened_hash, $flattened_names) = clone_features ($gene_array);
	my $result_string = "";
	foreach my $fullname (@$flattened_names) {
		my $strand = $flattened_hash->{$fullname}->{'strand'};
		my $start = $flattened_hash->{$fullname}->{'start'};
		my $end = $flattened_hash->{$fullname}->{'end'};
		my $regseq = $flattened_hash->{$fullname}->{'characters'};
		$result_string .= ">$fullname($strand)\t$start\t$end\n$regseq\n";
	}

	return $result_string;
}

sub write_sequin_tbl {
	my $gene_array = shift;
	my $name = shift;

	# print header
	my $result_string = ">Features\t$name\n";

	# start printing genes
	foreach my $gene (@$gene_array) {
		my $type = $gene->{'type'};
		foreach my $r (@{$gene->{'region'}}) {
			if ($r =~ /(\d+)\.\.(\d+)/) {
				# first, print the main feature:
				$result_string .= "$1\t$2\t$type\n";
				foreach my $q (keys %{$gene->{'qualifiers'}}) {
					my $qual = $gene->{qualifiers}->{$q};
					$result_string .= "\t\t\t$q\t$qual\n";
				}

				# then, print each feature contained.
				foreach my $feat (@{$gene->{'contains'}}) {
					$result_string .= Genbank::sequin_feature ($feat->{'region'}, $feat);
				}
			}
		}
	}
	return $result_string;
}

sub feature_table_from_genbank {
	my $gbfile = shift;
	my $gene_array = Genbank::simplify_genbank_array(Genbank::parse_genbank($gbfile));

	my @feature_array = ();
	foreach my $gene (@$gene_array) {
		if ($gene->{'type'} eq "gene") {
			# make an entry for the main gene
			my $feat_hash = {};
			$feat_hash->{'qualifiers'} = $gene->{'qualifiers'};
			$feat_hash->{'region'} = flatten_interval ($gene->{'region'});
			$feat_hash->{'type'} = $gene->{'type'};
			$feat_hash->{'contains'} = [];
			push @feature_array, $feat_hash;
			foreach my $feat (@{$gene->{'contains'}}) {
				# then make an entry for each subcomponent
				my $subfeat_hash = {};
				$subfeat_hash->{'qualifiers'} = $feat->{'qualifiers'};
				$subfeat_hash->{'region'} = $feat->{'region'};
				$subfeat_hash->{'type'} = $feat->{'type'};
				push @{$feat_hash->{'contains'}}, $subfeat_hash;
			}
		}
	}

	return \@feature_array;
}

sub sequence_for_interval {
	my $interval_str = shift;
	my $revcomp = shift;

	$interval_str =~ /(\d+)\.\.(\d+)/;
	my $start = $1;
	my $end = $2;
	my $geneseq = "";
	if ($start < $end) {
		(undef, $geneseq, undef) = Subfunctions::split_seq ($sequence, $start, $end);
	} else {
		(undef, $geneseq, undef) = Subfunctions::split_seq ($sequence, $end, $start);
		if ($revcomp == 1) {
			$geneseq = reverse_complement ($geneseq);
		}
	}
	return $geneseq;
}

sub sequin_feature {
	my $regions = shift;
	my $feature = shift;
	my $result_string = "";
	my $first_int = shift @$regions;
	$first_int =~ /(\d+)\.\.(\d+)/;
	$result_string = "$1\t$2\t$feature->{type}\n";
	foreach my $int (@$regions) {
		$int =~ /(\d+)\.\.(\d+)/;
		$result_string .= "$1\t$2\n";
	}
	foreach my $key (keys %{$feature->{'qualifiers'}}) {
		$result_string .= "\t\t\t$key\t$feature->{qualifiers}->{$key}\n";
	}
	return $result_string;
}

sub stringify_feature {
	my $regions = shift;
	my $feature = shift;
	my $result_string = "";

	my @features = ();
	foreach my $key (keys %{$feature->{'qualifiers'}}) {
		my $feat = "$key=$feature->{qualifiers}->{$key}";
		push @features, $feat;
	}
	$result_string = "$feature->{type}\t".join(",",@$regions)."\t".join("#",@features);
	return $result_string;
}

sub add_qualifiers_to_feature_hash {
	my $feature_hash = shift;
	my $feature_to_add = shift;

	foreach my $key (keys %{$feature_to_add->{'qualifiers'}}) {
		$feature_hash->{$key} = $feature_to_add->{'qualifiers'}->{$key};
	}
	return $feature_hash;
}

sub flatten_interval {
	my $int_array = shift;

	# calculate the largest extent of the main interval.
	my @locs = ();
	foreach my $r (@$int_array) {
		if ($r =~ /(\d+)\.\.(\d+)/) {
			push @locs, $1;
			push @locs, $2;
		}
	}
	# sort all the listed locations, make the start be the smallest possible val and the end be the largest possible val.
	@locs = sort {$a <=> $b} @locs;
	my $main_start = shift @locs;
	my $main_end = pop @locs;

	# was the original interval complemented?
	my $r = @$int_array[0];
	my ($loc1, $loc2) = split (/\.\./, $r);
	if ($loc1 < $loc2) {
		return "$main_start..$main_end";
	} else {
		return "$main_end..$main_start";
	}

}

sub parse_feature_desc {
	my $feat_desc_string = shift;
	my $gene_array_ptr = shift;

	my $feat_desc_hash = {};

	my ($feature_desc, @feature_quals) = split (/\+\//, $feat_desc_string);
	#parse feature_desc into name and interval.
	$feature_desc =~ s/\+//g;
	$feature_desc =~ /^(.{16})(.+)$/;
	my $type = $1;
	my $region = $2;
	$type =~ s/ *$//; # remove trailing blanks.

	$feat_desc_hash->{'type'} = $type;
	$feat_desc_hash->{'region'} = parse_interval($region);
	$feat_desc_hash->{'qualifiers'} = parse_qualifiers(\@feature_quals);
	if ($feat_desc_hash->{'type'} eq "gene") {
		$curr_gene = {};
		push @$gene_array_ptr, $curr_gene;
		$curr_gene->{'contains'} = ();
		$curr_gene->{'region'} = $feat_desc_hash->{'region'};
		$curr_gene->{'qualifiers'} = $feat_desc_hash->{'qualifiers'};
		$curr_gene->{'type'} = "gene";
	} else {
		if (within_interval($curr_gene->{'region'}, $feat_desc_hash->{'region'})) {
			# this feat_desc_hash belongs to the current gene.
			push @{$curr_gene->{'contains'}}, $feat_desc_hash;
		} else {
			$curr_gene = $feat_desc_hash;
			push @$gene_array_ptr, $feat_desc_hash;
		}
	}
	return $feat_desc_hash;
}

sub parse_interval {
	my $intervalstr = shift;
	my @regions = ();
	if ($intervalstr =~ /^complement\s*\((.+)\)/) {
		# this is a complementary strand feature.
		my $subregions = parse_interval($1);
		foreach my $subreg (@$subregions) {
			if ($subreg =~ /(\d+)\.\.(\d+)/) {
				unshift @regions, "$2..$1";
			}
		}
	} elsif ($intervalstr =~ /^join\s*\((.+)\)$/) {
		# this is a series of intervals
		my @subintervals = split(/,/, $1);
		foreach my $subint (@subintervals) {
			my $subregions = parse_interval($subint);
			push @regions, @$subregions;
		}
	} elsif ($intervalstr =~ /^order\s*\((.+)\)$/) {
		# this is a series of intervals, but there is no implication about joining them.
		my @subintervals = split(/,/, $1);
		my $max_interval = max_interval (\@subintervals);
		push @regions, @$max_interval;
	} elsif ($intervalstr =~ /(\d+)\.\.(\d+)/) {
		push @regions, "$intervalstr";
	}
	return \@regions;
}

sub parse_qualifiers {
	my $qualifiers_ref = shift;

	my @qualifiers = @$qualifiers_ref;
	my $feature_hash = {};
	while (@qualifiers > 0) {
		my $f = shift @qualifiers;
		$f =~ s/\+/ /g;
		if ($f =~ /(.+)=(.+)/) {
			my $key = $1;
			my $val = $2;
			if ($val =~ /"(.*)"/) {
				$val = $1;
			}
			if ($key eq "translation") {
				next;
				$val =~ s/ //g;
			}
			$feature_hash->{$key} = $val;
		} elsif ($f =~ /(.+)/) {
			my $key = $1;
			$feature_hash->{$key} = "";
		} else {
			print "haven't dealt with this: $f\n";
		}
	}
	return $feature_hash;
}

sub max_interval {
	my $regions = shift;
	# calculate the largest extent of the main interval.
	my @locs = ();
	my $strand = "+";
	foreach my $r (@$regions) {
		if ($r =~ /(\d+)\.\.(\d+)/) {
			push @locs, $1;
			push @locs, $2;
			if ($2 < $1) {
				$strand = "-";
			}
		}
	}
	# sort all the listed locations, make the start be the smallest possible val and the end be the largest possible val.
	@locs = sort {$a <=> $b} @locs;
	my $main_start = shift @locs;
	my $main_end = pop @locs;
	my @max_region = ();
	if ($strand eq "+") {
		push @max_region, "$main_start..$main_end";
	} else {
		push @max_region, "$main_end..$main_start";
	}
	return \@max_region;
}

sub within_interval {
	my $main_interval = shift;
	my $test_interval = shift;

	if (($main_interval eq "") || ($test_interval eq "")) {
		# if the interval we're testing for is blank, return 0.
		return 0;
	}

	if ((ref $test_interval) !~ /ARRAY/) {
		$test_interval = parse_interval($test_interval);
	}

	if ((ref $main_interval) !~ /ARRAY/) {
		$main_interval = parse_interval($main_interval);
	}

	# calculate the largest extent of the main interval.
	my @locs = ();
	foreach my $r (@$main_interval) {
		if ($r =~ /(\d+)\.\.(\d+)/) {
			push @locs, $1;
			push @locs, $2;
		}
	}
	# sort all the listed locations, make the start be the smallest possible val and the end be the largest possible val.
	@locs = sort {$a <=> $b} @locs;
	my $main_start = shift @locs;
	my $main_end = pop @locs;

	# do the same for the tested intervals.
	@locs = ();
	foreach my $r (@$test_interval) {
		if ($r =~ /(\d+)\.\.(\d+)/) {
			push @locs, $1;
			push @locs, $2;
		}
	}
	# sort all the listed locations, make the start be the smallest possible val and the end be the largest possible val.
	@locs = sort {$a <=> $b} @locs;
	my $test_start = shift @locs;
	my $test_end = pop @locs;
	if (($test_start >= $main_start) && ($test_end <= $main_end)) {
		return 1;
	}
	return 0;
}


return 1;

