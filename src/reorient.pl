#!/usr/bin/env perl
#===============================================================================
=pod

=head2

         FILE: reorient.pl

        USAGE: ./reorient.pl [options] infile.fasta

  DESCRIPTION: Try to "re-orient" (reverse+complement) sequences in
               infile.fasta based on blast hits.
               All sequences in input are blasted against the same
               reference sequence.
               A good strategy is to have as long and "representative"
               sequence as possible as reference.

      OPTIONS: -r,--ref=<ref.fas>   Use (first) sequence in file ref.fas as
                                    reference. If ref.fas is not
                                    given, the first sequence in input will
                                    be used as reference.
               -o,--out=<out.fas>   Write to out.fas instead of STDOUT.
               -f,--fail=<fail.fas> Sequences having 'No hit' with ref 
                                    sequence will be stored in fail.fas.
               -v,--version         Show version
               -h,--help            Show help
               -V,--verbose         Show info about output
               --noverbose          No verbose

 REQUIREMENTS: BioPerl, bl2seq (ncbi-blast+)

        NOTES: Current version uses blast "plus" and perl module
               Bio::Tools::Run::StandAloneBlastPlus.
               Warnings from bl2seq may be printed to stderr.

       AUTHOR: Johan Nylander

      COMPANY: NRM

      VERSION: 1.0.1

      CREATED: 10/07/2009 01:04:38 PM CEST

     REVISION: tis 16 jan 2024 11:02:31

      LICENSE: Copyright (c) 2009-2024 Johan Nylander

               Permission is hereby granted, free of charge, to any person
               obtaining a copy of this software and associated documentation
               files (the "Software"), to deal in the Software without
               restriction, including without limitation the rights to use,
               copy, modify, merge, publish, distribute, sublicense, and/or
               sell copies of the Software, and to permit persons to whom the
               Software is furnished to do so, subject to the following
               conditions:

               The above copyright notice and this permission notice shall be
               included in all copies or substantial portions of the Software.

               THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
               EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
               OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
               NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
               HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
               WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
               FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
               OTHER DEALINGS IN THE SOFTWARE.

=cut
#===============================================================================

use strict;
use warnings;

use Data::Dumper;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlastPlus;
use Getopt::Long;
Getopt::Long::Configure("no_ignore_case", "no_auto_abbrev");

## File globals
my $version        = '1.0.1';
my $ref_file       = q{};
my $out_file       = q{};
my $output_fh      = q{};
my $fail_file      = q{};
my $fail_fh        = q{};
my $blast_out_file = q{};
my $format         = 'fasta';
my $did_rc         = 0;
my $VERBOSE        = 1;

exec("perldoc", $0) unless (@ARGV);

## Options
my $res = GetOptions(
    'r|ref=s'    => \$ref_file,
    'o|out=s'    => \$out_file,
    'f|fail=s'   => \$fail_file,
    'V|verbose!' => \$VERBOSE,
    'v|version'  => sub { print "$version\n"; exit(0); },
    'h|help'     => sub {exec("perldoc", $0); exit(0);}
);

if ($out_file) {
    $output_fh = Bio::SeqIO->newFh(-format => $format, -file => ">$out_file");
}
else {
    $output_fh = Bio::SeqIO->newFh(-format => $format);
}

if ($fail_file) {
    $fail_fh = Bio::SeqIO->new(-format => $format, -file => ">$fail_file");
}

## Read from standard input or the input filenames
my $stream = Bio::SeqIO->newFh(-format => $format, -fh => \*ARGV);

## Open ref seq from file, otherwise use the first seq in input
my $ref_seq_obj = '';
if ($ref_file) {
    my $ref_file_io_object = Bio::SeqIO->new(-file => $ref_file, -format => $format);
    $ref_seq_obj = $ref_file_io_object->next_seq();
}
else {
    $ref_seq_obj = <$stream>; # capture first seq in file
    $output_fh->print($ref_seq_obj); # if no ref file is used, print also the first seq
}

## Iterate over the rest of the sequences in stream
while (my $query_seq_obj = <$stream>) {

    my $ref_input = Bio::Seq->new(-id => "reference", -seq => $ref_seq_obj->seq);
    my $query_input = Bio::Seq->new(-id => "sequence", -seq => $query_seq_obj->seq);

    my $factory = Bio::Tools::Run::StandAloneBlastPlus->new();

    my $blast_result = $factory->bl2seq(
        -method => 'blastn',
        -query => $query_input,
        -subject => $ref_input
    );

    if (my $hit_obj = $blast_result->next_hit) {

        ## there may be >1 hsp, but I'm only looking at the first
        my $hsp_obj = $hit_obj->next_hsp;

        ## Reverse+complement if hit strand is '-1'
        if ($hsp_obj->strand('hit') < 1) {
            if ($query_seq_obj->alphabet eq 'dna') {
                my $revcomp_query_seq_obj = $query_seq_obj->revcom();
                $query_seq_obj->seq($revcomp_query_seq_obj->seq);
                my $id = $query_seq_obj->display_id();
                $id = "$id\@revcomp\@";
                $query_seq_obj->display_id($id);
                $output_fh->print($query_seq_obj);
                $did_rc = 1;
            }
        }
        else {
            $output_fh->print($query_seq_obj);
        }
    }
    else {
        my $fail_id = $query_seq_obj->display_id();
        if ($VERBOSE) {
            print STDERR "Warning: No hits for $fail_id\n";
        }
        if ($fail_file) {
            $fail_fh->write_seq($query_seq_obj);
        }
    }
    $factory->cleanup();
}

if ($did_rc) {
    if ($VERBOSE) {
        print STDERR "Note: At least one sequence have been reverse-complemented. \'\@revcomp@\' have been appended to the fasta header\n";
    }
}

if ($fail_file) {
    if (-z $fail_file) {
        if ($VERBOSE) {
            print STDERR "Note: All sequences had hits with reference. Removing fail file \"$fail_file\"\n";
        }
        unlink $fail_file;
    }
    else {
        if ($VERBOSE) {
            print STDERR "Warning: Some sequences did not have hits with reference. See fail file \"$fail_file\"\n";
        }
    }
}

exit(0);

