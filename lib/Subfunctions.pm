package Subfunctions;
use strict;
use File::Temp qw(tempfile);
use FindBin;
use lib "$FindBin::Bin/lib";
use Blast;
use Genbank;
use Data::Dumper;

BEGIN {
	require Exporter;
	# set the version for version checking
	our $VERSION     = 1.00;
	# Inherit from Exporter to export functions and variables
	our @ISA         = qw(Exporter);
	# Functions and variables which are exported by default
	our @EXPORT      = qw();
	# Functions and variables which can be optionally exported
	our @EXPORT_OK   = qw(parse_fasta blast_to_genbank align_regions_to_reference align_hits_to_ref split_seq reverse_complement);
}

my $debug = 0;

sub set_debug {
	$debug = shift;
}

sub debug {
	my $msg = shift;
	if ($debug) {
		print $msg;
	}
}

=head1

B<(\%taxa, \@taxanames) parse_fasta ( Filehandle $inputfile )>

Given a fasta file as input, returns a hash containing all the sequences, keyed by the
values of the taxanames array.

$inputfile:   fasta file to parse.

=cut


sub parse_fasta {
	my $fastafile = shift;
	my $no_substitution = shift;

	my $taxa = {};
	my @taxanames = ();
	my $length = 0;
	open fileIN, "<:crlf", "$fastafile" or die "couldn't parse fasta file $fastafile";

	my $input = readline fileIN;
	my $taxonlabel = "";
	my $sequence = "";
	while (defined $input) {
		if ($input =~ /^>(.+)$/) {
			$taxonlabel = $1;
			$taxonlabel =~ s/\s+$//;
			if (!defined $no_substitution) {
				$taxonlabel =~ s/[\s\\\/,;\-]+/_/g;
			}
			if (exists $taxa->{$taxonlabel}) {
				$taxa->{$taxonlabel} = "";
			} else {
				push @taxanames, $taxonlabel;
			}
			if ($length > 0) {
				# we are at the next taxon; push the last one onto the taxon array.
				$taxa->{'length'} = $length;
				$length = 0;
			}
		} else {
			if ($input =~ /^\s*(.+)\s*$/) {
				$taxa->{$taxonlabel} .= $1;
				$length += length($1);
			}
		}
		$input = readline fileIN;
	}

	close (fileIN);
	return $taxa, \@taxanames;
}

# returns an array of gene hashes:
# gene->{name}
# gene->{strand}
# gene->{type}
# gene->{contains} = an array of regions
#	[ (start, end), (start, end) ]

sub blast_to_genbank {
	my $gbfile = shift;
	my $fastafile = shift;

	my $gene_array = Genbank::simplify_genbank_array(Genbank::parse_genbank($gbfile));
	my $refseq = Genbank::get_sequence();
	my ($ref_hash, $region_array) = Genbank::clone_features($gene_array);
	my ($query_hash, $query_array) = parse_fasta($fastafile);
	my $queryseq = $query_hash->{@$query_array[0]};

	# look for regions too small to blast accurately:
	my $tiny_regions = {};
	my $tiny_region_extension_length = 20;
	foreach my $region (@$region_array) {
		my $start = $ref_hash->{$region}->{'start'};
		my $end = $ref_hash->{$region}->{'end'};
		if ($end - $start < 10) {
			$tiny_regions->{$region}->{'characters'} = $ref_hash->{$region}->{'characters'};
			$start -= $tiny_region_extension_length;
			$end += $tiny_region_extension_length;
			my (undef, $seq, undef) = Subfunctions::split_seq ($refseq, $start, $end);
			$ref_hash->{$region}->{'characters'} = $seq;
		}
	}

 	my ($fastafh, $subjectfasta) = tempfile();
# 	my $subjectfasta = "temp.fasta";
	open my $fastafh, ">", $subjectfasta;

	my @main_gene_array = ();
	foreach my $ref (@$region_array) {
		# if 'features' exists, it's a main gene.
		if (exists $ref_hash->{$ref}->{'qualifiers'}) {
			push @main_gene_array, $ref;
			# if it is a main gene, only push the sequence into the fasta file if it doesn't have subcomponents.
			if (@{$ref_hash->{$ref}->{'contains'}} == 0) {
				print $fastafh ">$ref\t$ref_hash->{$ref}->{'start'}\t$ref_hash->{$ref}->{'end'}\n$ref_hash->{$ref}->{'characters'}\n";
			}
		} else {
			print $fastafh ">$ref\t$ref_hash->{$ref}->{'start'}\t$ref_hash->{$ref}->{'end'}\n$ref_hash->{$ref}->{'characters'}\n";
		}
	}
	close $fastafh;

 	my (undef, $blastfile) = tempfile();
# 	my $blastfile = "temp.xml";
	system("blastn -query $fastafile -subject $subjectfasta -outfmt 5 -out $blastfile -word_size 10");

	# choose the best hits:
	my $hit_array = Blast::parse_xml ($blastfile);
	foreach my $hit (@$hit_array) {
		my $subj = $hit->{'subject'}->{'name'};

		# keep the best-scoring hits.
		my @hsps = sort Blast::sort_hsps_by_score @{$hit->{'hsps'}};
		my @best_hsps = ();
		foreach my $hsp (@hsps) {
			if ($hsp->{'bit-score'} > 30) {
				push @best_hsps, $hsp;
			} else {
				last;
			}
		}
		my $best_hit = $best_hsps[0];

		# fix the lengths back for tiny regions.
		if (exists $tiny_regions->{$subj}) {
			$best_hit->{'hit-from'} += $tiny_region_extension_length;
			$best_hit->{'hit-to'} -= $tiny_region_extension_length;
			$best_hit->{'query-from'} += $tiny_region_extension_length;
			$best_hit->{'query-to'} -= $tiny_region_extension_length;
		}

		$best_hit->{'align-len'} = $best_hit->{'hit-to'} - $best_hit->{'hit-from'} + 1;
		$ref_hash->{$subj}->{'align-len'} = ($ref_hash->{$subj}->{'end'} - $ref_hash->{$subj}->{'start'} + 1);

		my (undef, $refgeneseq, undef) = Subfunctions::split_seq ($refseq, $ref_hash->{$subj}->{'start'}, $ref_hash->{$subj}->{'end'});
		$ref_hash->{$subj}->{'reference'} = $refgeneseq;
		$ref_hash->{$subj}->{'refstart'} = delete $ref_hash->{$subj}->{'start'};
		$ref_hash->{$subj}->{'refend'} = delete $ref_hash->{$subj}->{'end'};
		$ref_hash->{$subj}->{'complete'} = 0;
		$ref_hash->{$subj}->{'gaps'} = "none";
		$ref_hash->{$subj}->{'hsps'} = \@best_hsps;
	}

	# check to see if each of the main genes is actually complete:
	foreach my $main_gene (@main_gene_array) {
		$ref_hash->{$main_gene}->{'complete'} = 1;

		# if there weren't any components, check to see if there were hits for the main gene
		if (@{$ref_hash->{$main_gene}->{'contains'}} == 0) {
			if (@{$ref_hash->{$main_gene}->{'hsps'}} == 0) {
				$ref_hash->{$main_gene}->{'complete'} = 0;
			} else {
				my $best_hit = @{$ref_hash->{$main_gene}->{'hsps'}}[0];
				# copy the best hit values back into the reference array:
				my @best_one = ($best_hit);
				$ref_hash->{$main_gene}->{'start'} = $best_hit->{'query-from'};
				$ref_hash->{$main_gene}->{'end'} = $best_hit->{'query-to'};
				if ($best_hit->{'align-len'} == $ref_hash->{$main_gene}->{'align-len'}) {

					my $hseq_gaps = $best_hit->{'hseq'};
					my $qseq_gaps = $best_hit->{'qseq'};

					$ref_hash->{$main_gene}->{'gaps'} = ($hseq_gaps =~ tr/[AGCT]//dr) . " " . ($qseq_gaps =~ tr/[AGCT]//dr);
					if ((($hseq_gaps =~ /-/) ne "1") && (($qseq_gaps =~ /-/) ne "1")) {
						$ref_hash->{$main_gene}->{'complete'} = 1;
					}
				}
			}
		} else {
			foreach my $subj (@{$ref_hash->{$main_gene}->{'contains'}}) {
				my $best_hit = @{$ref_hash->{$subj}->{'hsps'}}[0];
				# copy the best hit values back into the reference array:
				my @best_one = ($best_hit);
				$ref_hash->{$subj}->{'hsps'} = \@best_one;
				$ref_hash->{$subj}->{'start'} = $best_hit->{'query-from'};
				$ref_hash->{$subj}->{'end'} = $best_hit->{'query-to'};
				if ($best_hit->{'align-len'} == $ref_hash->{$subj}->{'align-len'}) {

					# need to check for whether or not the gaps are in threes.
					my $hseq_gaps = $best_hit->{'hseq'};
					my $qseq_gaps = $best_hit->{'qseq'};
					$hseq_gaps =~ s/---//g;
					$qseq_gaps =~ s/---//g;

					$ref_hash->{$subj}->{'gaps'} = ($hseq_gaps =~ tr/[AGCT]//dr) . " " . ($qseq_gaps =~ tr/[AGCT]//dr);
					if ((($hseq_gaps =~ /-/) ne "1") && (($qseq_gaps =~ /-/) ne "1")) {
						$ref_hash->{$subj}->{'complete'} = 1;
					}
				}

				# if one component is not complete, the whole gene is not complete.
				if ($ref_hash->{$subj}->{'complete'} == 0) {
					$ref_hash->{$main_gene}->{'complete'} = 0;
				}
			}
		}
	}
	return ($ref_hash, \@main_gene_array);
}

=head1

B<align_regions_to_reference>

Takes in a hash and index array from blast_to_genbank or similar as well as a Genbank record.
Reindexes the positions of the features in the Genbank record to the hashed values.
Returns an array of hashes for each feature:
    hash->{'type'} = the feature type
    hash->{'qualifiers'} = Genbank qualifiers, in a hash
    hash->{'region'} = an array of intervals (as strings of xx..xx format)
    hash->{'contains'} = the subfeatures.

=cut

sub align_regions_to_reference {
	my $region_hash = shift;
	my $region_array = shift;
	my $refgbfile = shift;

	my $ref_gene_array = Genbank::feature_table_from_genbank ($refgbfile);

	# clone the genes from ref_gene_array, replace with the info from the matched_gene_array
	my @final_gene_array = ();
	my $new_gene_count = 0;
	my $dest_gene_count = 0;
	my $new_gene = $region_hash->{@$region_array[$new_gene_count++]};
	my $dest_gene = @$ref_gene_array[$dest_gene_count++];

	while (defined $dest_gene) {
			my $destgenename = "";
			if ($dest_gene->{'qualifiers'}->{'locus_tag'}) {
				$destgenename = $dest_gene->{'qualifiers'}->{'locus_tag'};
			} else {
				$destgenename = $dest_gene->{'qualifiers'}->{'gene'};
			}
			my $newgenename = "";
			if ($new_gene->{'qualifiers'}->{'locus_tag'}) {
				$newgenename = $new_gene->{'qualifiers'}->{'locus_tag'};
			} else {
				$newgenename = $new_gene->{'qualifiers'}->{'gene'};
			}


		if ($newgenename ne $destgenename) {
			$dest_gene = @$ref_gene_array[$dest_gene_count++];
			next;
		}

		my $final_gene = {};
		# dest gene has:
		#	type:
		$final_gene->{'type'} = $dest_gene->{'type'};

		#	qualifiers:
		foreach my $q (keys %{$dest_gene->{'qualifiers'}}) {
			$final_gene->{'qualifiers'}->{$q} = $dest_gene->{'qualifiers'}->{$q};
		}

		# dest gene has a contains array: each of these has a subcomponent matching a contains piece in new gene
		$final_gene->{'contains'} = ();
		my @final_regions = ();

		foreach my $subcomponent (@{$dest_gene->{'contains'}}) {
			# each subcomponent will have qualifiers, type, region array
			my $newcontains = {};
			push @{$final_gene->{'contains'}}, $newcontains;
			$newcontains->{'qualifiers'} = $subcomponent->{'qualifiers'};
			$newcontains->{'type'} = $subcomponent->{'type'};
			$newcontains->{'region'} = \@final_regions;
			foreach my $subcomp (@{$new_gene->{'contains'}}) {
				my $reg = "$region_hash->{$subcomp}->{'start'}..$region_hash->{$subcomp}->{'end'}";
				if ($region_hash->{$subcomp}->{'strand'} eq "-") {
					$reg = "$region_hash->{$subcomp}->{'end'}..$region_hash->{$subcomp}->{'start'}";
				} 
				
				push @final_regions, $reg;
			}
		}

		# if there weren't any regions in subcomponents, make the region out of the main gene region
		if (@final_regions == 0) {
			my $reg = "$new_gene->{'start'}..$new_gene->{'end'}";
			if ($new_gene->{'strand'} eq "-") {
				$reg = "$new_gene->{'end'}..$new_gene->{'start'}";
			} 
			
			push @final_regions, $reg;
		}
		#	region: this should be a stringified max interval from the new gene
		$final_gene->{'region'} = Genbank::max_interval(\@final_regions);

		push @final_gene_array, $final_gene;
		$new_gene = $region_hash->{@$region_array[$new_gene_count++]};

# 		if (!defined $new_gene) {
# 			last;
# 		}
	}
	return \@final_gene_array;
}

sub align_hits_to_ref {
	my $ref_hash = shift;
	my $gene_name = shift;

	my @bits = ();
	if (@{$ref_hash->{$gene_name}->{'contains'}} > 0) {
		foreach my $bit (@{$ref_hash->{$gene_name}->{'contains'}}) {
			push @bits, $bit;
		}
	} else {
		push @bits, $gene_name;
	}
	my $result_string = "";
	foreach my $bit (@bits) {
		$result_string = "  REFERENCE $bit:\n";
		my $hit_hash = $ref_hash->{$bit};
		my @hits = sort Blast::sort_hsps_by_hit_start @{$hit_hash->{'hsps'}};
		$result_string .= " " x (11 - (length $hit_hash->{'refstart'})) . "$hit_hash->{'refstart'} ";
		$result_string .= "$hit_hash->{'reference'} $hit_hash->{'refend'}\n";
		$result_string .= "  MATCHES:\n";
		my $offset = $hit_hash->{'refstart'} - 1;
		foreach my $hit (@hits) {
			if ($hit->{'hit-from'} > $hit->{'hit-to'}) {
				next;
			}
			my $front_pad = ($hit->{'hit-from'} - 1);
			my $back_pad = (length($hit_hash->{'reference'}) - $hit->{'hit-to'});

			my $hseq = '-'x $front_pad . $hit->{'hseq'} . '-' x $back_pad;
			$result_string .= "Ref" . " " x (8 - (length ($hit->{'hit-from'} + $offset))) . ($hit->{'hit-from'} + $offset) . " ";
			$result_string .= "$hseq " . ($hit->{'hit-to'} + $offset) . "\n";

			my $midline = '-'x $front_pad . $hit->{'midline'} . '-' x $back_pad;
			$result_string .= " " x 12;
			$result_string .= "$midline\n";

			my $qseq = '-'x $front_pad . $hit->{'qseq'} . '-' x $back_pad;
			$result_string .= "New" . " " x (8 - (length $hit->{'query-from'})) . "$hit->{'query-from'} ";
			$result_string .= "$qseq $hit->{'query-to'}\n\n";
		}
	}
 	return $result_string;
}

=head1

B<String ($startseq, $regionseq, $endseq) split_seq ( String $seq, int $start, int $end )>

=cut
### taken from aTRAM:

sub split_seq {
    my $seq = shift;
    my $start = shift;
    my $end = shift;
    my $max = 30000;
	my $seqlen = length ($seq);
	my $startseq = "";
	my $regionseq = "";
	my $endseq = "";

	my $currstart = $start-1;
	my $currend = $end;
	while ($currstart > $max) {
			$seq =~ /^(.{$max})(.*)$/;
			$startseq .= $1;
			$seq = $2;
			$currstart -= $max;
			$currend -= $max;
	}
	if ($currstart > 0) {
			$seq =~ /^(.{$currstart})(.*)$/;
			$startseq .= $1;
			$seq = $2;
			$currstart = 1;
			$currend = $end - (length ($startseq));
	}

	my $regionsize = $end - $start + 1;
	while ($regionsize > $max) {
			$seq =~ /^(.{$max})(.*)$/;
			$regionseq .= $1;
			$seq = $2;
			$currstart -= $max;
			$currend -= $max;
			$regionsize -= $max;
	}
	if ($regionsize > 0) {
			$seq =~ /^(.{$regionsize})(.*)$/;
			$regionseq .= $1;
			$endseq = $2;
	}
	return ($startseq, $regionseq, $endseq);
}

=head1

B<String reverse_complement ( String $charstr )>

Convenience function to return the reverse complement of a sequence.

$charstr:   sequence to revcomp.

=cut


sub reverse_complement {
	my $charstr = shift;

	# reverse the DNA sequence
	my $revcomp = reverse($charstr);

	# complement the reversed DNA sequence
	$revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
	return $revcomp;
}

# must return 1 for the file overall.
1;
