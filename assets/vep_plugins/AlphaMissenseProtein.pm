=head1 LICENSE

Copyright the lrsomatic authors.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

 https://github.com/IntGenomicsLab/lrsomatic

=cut

=head1 NAME

 AlphaMissenseProtein

=head1 SYNOPSIS

 # default columns: am_pathogenicity and am_class
 ./vep -i variants.vcf --dir_plugins . \
   --plugin AlphaMissenseProtein,file=/path/to/table.tsv.gz

 # pick columns, or cols=all. VEP splits plugin parameters on commas, so the
 # column list is '&' separated
 ./vep -i variants.vcf --dir_plugins . \
   --plugin AlphaMissenseProtein,file=/path/to/table.tsv.gz,cols=am_class&uniprot_acc

=head1 DESCRIPTION

 Annotates amino-acid-changing variants with AlphaMissense pathogenicity values
 looked up in PROTEIN space rather than genomic space, which is what makes them
 reachable on T2T-CHM13. AlphaMissense is published in GRCh37/GRCh38 coordinates
 only, but the score is a function of (protein sequence, amino-acid
 substitution), so keyed on gene symbol and substitution it is
 assembly-independent and needs no liftover.

 The table is that release re-keyed from UniProt accession onto gene symbol,
 indexed with `tabix -s 1 -b 2 -e 2` -- tabix's "sequence" column is just a
 string, which is what lets a protein-space table be range-queried at all. See
 CITATIONS.md for how it was derived.

 A row is used only when both its reference and its alternate amino acid equal
 what VEP computed for the transcript at hand, so where the CHM13 protein
 differs from the GRCh38 one it returns no score rather than a score for the
 wrong substitution. AlphaMissenseProtein_match reports which happened:
   gene_aa      - matched on gene symbol and amino-acid substitution
   aa_mismatch  - the gene and position exist but no row has this substitution
   not_found    - the gene symbol and position are absent from the table
   no_gene      - VEP produced no gene symbol for this transcript

 Requires --symbol (implied by --everything), since the join key is the gene
 symbol.

 ATTRIBUTION: the data is the AlphaMissense Database, Copyright (2023) DeepMind
 Technologies Limited, licensed CC BY 4.0. CHANGES WERE MADE: re-keyed from
 UniProt accession to gene symbol and reshaped into a tabix-indexed table, with
 no value altered. Provided for theoretical modelling only, not as a substitute
 for professional medical advice. Cite https://doi.org/10.1126/science.adg7492.

=cut

package AlphaMissenseProtein;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

# Key columns of the lookup table. Never emitted as annotation -- they are the
# join key and the guard, not results.
my @KEY_COLS = qw(gene aapos aaref aaalt);

my @DEFAULT_COLS = qw(am_pathogenicity am_class);

my %DESCRIPTIONS = (
  am_pathogenicity => 'AlphaMissense pathogenicity (protein-space lookup), 0-1',
  am_class         => 'AlphaMissense class: likely_pathogenic (>0.564), likely_benign (<0.34), ambiguous',
  uniprot_acc      => 'UniProt accession the AlphaMissense numbering came from',
);

sub new {
  my $class = shift;
  my $self = $class->SUPER::new(@_);

  # No region padding: we want exact amino-acid positions, and disabling the
  # region cache means get_data() returns unfiltered rows we filter ourselves.
  $self->expand_left(0);
  $self->expand_right(0);

  my $param_hash = $self->params_to_hash();

  # Accept file= for consistency with the stock AlphaMissense plugin, but also
  # a bare positional path.
  my $file = $param_hash->{file};
  if (!$file) {
    my @positional = grep { !/=/ } @{$self->params};
    $file = shift @positional;
  }

  die "ERROR: AlphaMissenseProtein needs a path to the protein score table, e.g. --plugin AlphaMissenseProtein,file=/path/to/table.tsv.gz\n"
    unless $file;

  $self->add_file($file);

  $self->{headers} = $self->_read_table_header($file);
  my %is_header = map { $_ => 1 } @{$self->{headers}};
  my %is_key    = map { $_ => 1 } @KEY_COLS;

  for my $k (@KEY_COLS) {
    die "ERROR: protein score table $file is missing the required key column '$k'\n"
      unless $is_header{$k};
  }

  my @want;
  if ($param_hash->{cols} && $param_hash->{cols} eq 'all') {
    @want = grep { !$is_key{$_} } @{$self->{headers}};
  }
  elsif ($param_hash->{cols}) {
    my @requested = split /&/, $param_hash->{cols};
    my @invalid = grep { !$is_header{$_} } @requested;
    warn "WARNING: AlphaMissenseProtein: columns not in table header, ignored: " .
         join(',', @invalid) . "\n" if @invalid;
    @want = grep { $is_header{$_} && !$is_key{$_} } @requested;
  }
  else {
    @want = grep { $is_header{$_} } @DEFAULT_COLS;
  }

  die "ERROR: AlphaMissenseProtein: no valid columns selected. Table has:\n" .
      join(',', @{$self->{headers}}) . "\n" unless @want;

  $self->{cols} = \@want;

  return $self;
}

sub feature_types { return ['Transcript']; }

sub get_header_info {
  my $self = shift;
  my %info = map {
    $_ => ($DESCRIPTIONS{$_} || "$_ from the AlphaMissense protein-space table")
  } @{$self->{cols}};
  $info{AlphaMissenseProtein_match} =
    'How the protein-space lookup resolved: gene_aa (matched), aa_mismatch, not_found, or no_gene. ' .
    'Data: AlphaMissense Database (c) 2023 DeepMind Technologies Limited, CC BY 4.0; ' .
    're-keyed from UniProt accession to gene symbol, values unaltered';
  return \%info;
}

# The table's header is its first line, written with a leading '#' and indexed
# with `tabix -c '#'`. Prefer `tabix -H`; fall back to reading the bgzf
# directly, since bgzip output is gzip-compatible.
sub _read_table_header {
  my ($self, $file) = @_;

  for my $cmd ("tabix -H $file 2>/dev/null |", "gzip -dc $file 2>/dev/null |") {
    my @fields;
    if (open my $fh, $cmd) {
      while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next unless $line =~ /^#/;
        $line =~ s/^#//;
        @fields = split /\t/, $line;
        last;
      }
      close $fh;
    }
    return \@fields if @fields;
  }

  die "ERROR: could not read a '#'-prefixed header line from $file\n";
}

sub parse_data {
  my ($self, $line) = @_;

  $line =~ s/\r$//;
  return undef if $line =~ /^#/;

  my @split = split /\t/, $line;
  my %data;
  my $headers = $self->{headers};
  $data{$headers->[$_]} = $split[$_] for 0 .. $#$headers;

  return \%data;
}

sub get_start { return $_[1]->{aapos}; }
sub get_end   { return $_[1]->{aapos}; }

sub run {
  my ($self, $tva) = @_;

  my $tv = $tva->transcript_variation;
  return {} unless $tv;

  # A single-residue amino-acid substitution is the only thing this table keys on.
  my $pos = $tv->translation_start;
  my $end = $tv->translation_end;
  return {} unless defined $pos && defined $end && $pos == $end;

  my $pep = $tva->pep_allele_string;
  return {} unless defined $pep;

  my ($ref_aa, $alt_aa) = split m{/}, $pep;
  # No slash means synonymous: nothing to look up.
  return {} unless defined $ref_aa && defined $alt_aa && $ref_aa ne $alt_aa;

  my $tr = $tva->transcript;
  my $symbol = $tr->{_gene_symbol} || $tr->{_gene_hgnc};
  return { AlphaMissenseProtein_match => 'no_gene' }
    unless defined $symbol && $symbol ne '';

  # Query a small window and filter exactly. get_data() is 1-based on the
  # Bio::DB::HTS path and 0-based on the Tabix.pm path, so padding then
  # exact-matching aapos is correct on either backend; it also avoids passing a
  # start of 0 for the initiator methionine, which get_data() rejects.
  my $qs = $pos > 2 ? $pos - 2 : 1;
  my $rows = $self->get_data($symbol, $qs, $pos + 1);

  return { AlphaMissenseProtein_match => 'not_found' } unless $rows && @$rows;

  my $hit;
  for my $row (@$rows) {
    next unless defined $row->{aapos} && $row->{aapos} eq "$pos";
    next unless defined $row->{aaref} && $row->{aaref} eq $ref_aa;
    next unless defined $row->{aaalt} && $row->{aaalt} eq $alt_aa;
    $hit = $row;
    last;
  }

  # The position exists in the table but not this substitution: the local
  # protein and the protein AlphaMissense was numbered against disagree here, so
  # we deliberately return no score.
  return { AlphaMissenseProtein_match => 'aa_mismatch' } unless $hit;

  my %return = (AlphaMissenseProtein_match => 'gene_aa');
  for my $col (@{$self->{cols}}) {
    my $v = $hit->{$col};
    next unless defined $v && $v ne '.' && $v ne '';
    $return{$col} = $v;
  }

  return \%return;
}

1;
