#!/usr/bin/env perl

use Getopt::Long;
use Pod::Usage;
use File::Temp qw (tempfile);
use FindBin;
use lib "$FindBin::Bin/lib";
use Blast;
use Genbank;
use Subfunctions qw (parse_fasta blast_to_genbank align_regions_to_reference align_hits_to_ref);
use Data::Dumper;

my $help = 0;
my $outfile = "output";
my $gbfile = "";
my $fastafile = "";
my $orgname = "";
my $samplename = "";

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

GetOptions ('reference|gb|genbank=s' => \$gbfile,
			'fastafile=s' => \$fastafile,
			'outfile=s' => \$outfile,
			'organism=s' => \$orgname,
			'sample=s' => \$samplename,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 2);
}

if ($gbfile !~ /\.gbf*$/) {
	print "reference file needs to be a fully annotated Genbank file.\n";
	exit;
}

if ($gbfile eq "") {
	print "need to supply a Genbank reference file (-reference).\n";
	exit;
}

if ($fastafile eq "") {
	print "need to supply a fasta-formatted plastome sequence to annotate (-fasta).\n";
	exit;
}

my ($fastahash, $fastaarray) = parse_fasta($fastafile);

# there should be only one key, so just one name.
if ($samplename eq "") {
	$samplename = @$fastaarray[0];
}
my $queryseq = $fastahash->{@$fastaarray[0]};

print "Comparing reference genbank $gbfile to new plastome $fastafile\n";
my ($result_hash, $main_gene_array) = blast_to_genbank ($gbfile, $fastafile);
my @finished_array = ();
my $missing_results = "";
my $num_missing = 0;
print "Checking for pseudogenes and missing hits\n";
foreach my $main_gene (@$main_gene_array) {
	my $gene_features = $result_hash->{$main_gene}->{'qualifiers'};
	if ($result_hash->{$main_gene}->{'complete'} == 0) {
		print "MISSING $main_gene\n";
		$missing_results .= "MISSING $main_gene " . $result_hash->{$main_gene}->{'strand'} . " " . $result_hash->{$main_gene}->{'gaps'} . "\n" . align_hits_to_ref ($result_hash, $main_gene);
		$num_missing++;
	} else {
		push @finished_array, $main_gene;
		next;
	}


}

my $genbank_header = "$samplename [organism=$orgname][moltype=Genomic DNA][location=chloroplast][topology=Circular][gcode=11]";
open FASTA_FH, ">", "$outfile.fsa";
print FASTA_FH ">$genbank_header\n$queryseq\n";
close FASTA_FH;

print "Aligning matched regions to reference genbank file\n";
my $gene_array = align_regions_to_reference ($result_hash, \@finished_array, $gbfile);
# need to annotate inverted repeats
my ($fh, $refblast) = tempfile();
system("blastn -query $fastafile -subject $fastafile -outfmt 5 -out $refblast -evalue 1e-200");

my $self_array = Blast::parse_xml ("$refblast");
my @irs = ();
foreach my $hit (@$self_array) {
	my @hsps = sort Blast::sort_regions_by_start @{$hit->{'hsps'}};
	foreach my $hsp (@hsps) {
		my $querylen = $hsp->{'query-to'} - $hsp->{'query-from'};
		# IRs are between 10,000 and 50,000 bp and are inverted.
		if (($querylen < 50000) && ($querylen > 10000) && ($hsp->{'hit-frame'} == -1)) {
			push @irs, $hsp;
		}
	}
}

if (@irs > 2) {
	print "Warning! There seem to be more than two inverted repeats (".@irs." found). Are you sure this is a plastome sequence?\n";
}
@irs = sort Blast::sort_hsps_by_query_start @irs;

# all of the IR information is in one of the IRs, so shift.
my $ir = shift @irs;

my $ira = $ir->{'hit-to'} . ".." . $ir->{'hit-from'};
my $ira_hash = {};
$ira_hash->{'type'} = "repeat_region";
$ira_hash->{'qualifiers'}->{'note'} = "inverted repeat A";
$ira_hash->{'qualifiers'}->{'rpt_type'} = "inverted";
$ira_hash->{'region'} = ();
push @{$ira_hash->{'region'}}, $ira;
$ira_hash->{'contains'} = ();
push @$gene_array, $ira_hash;

my $irb = $ir->{'query-from'} . ".." . $ir->{'query-to'};
my $irb_hash = {};
$irb_hash->{'type'} = "repeat_region";
$irb_hash->{'qualifiers'}->{'note'} = "inverted repeat B";
$irb_hash->{'qualifiers'}->{'rpt_type'} = "inverted";
$irb_hash->{'region'} = ();
push @{$irb_hash->{'region'}}, $irb;
$irb_hash->{'contains'} = ();
push @$gene_array, $irb_hash;

open TBL_FH, ">", "$outfile.tbl";
print TBL_FH Genbank::write_sequin_tbl ($gene_array, $genbank_header);
close TBL_FH;

open MISSING_FH, ">", "$outfile.results.txt";
print MISSING_FH "Aligned " . @finished_array . " genes\n";
print MISSING_FH "Missing $num_missing genes\n\n" . $missing_results;

close MISSING_FH;


__END__

=head1 NAME

plann.pl

=head1 SYNOPSIS

plann.pl -reference gbfile.gb -fasta denovoplastome.fasta -out outfile [-organism "Genus species"] [-sample samplename]

=head1 OPTIONS

  -reference:       a well-annotated plastome reference sequence, in genbank file format
  -fastafile:       the plastome sequence to be annotated, in fasta format
  -outfile:         the output name (default is "output")
  -organism:        [optional: scientific name for Genbank annotation]
  -sample:          the name of the plastome sample (default is the name in the fasta file)

=head1 DESCRIPTION

Plann uses a reference plastome Genbank file as a template to annotate a related plastome sequence.

=cut
