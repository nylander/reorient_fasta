# Re-orient fasta sequences

- Last modified: tis jan 16, 2024  06:41
- Sign: nylander

## Description

Try to correct the direction of DNA sequences in fasta (by reverse complement)
in relation to a reference sequence.

More specifically, try to "re-orient" (reverse+complement) sequences based on
their blast hits (using bl2seq, blastn). All sequences in input are blasted
against the same reference sequence. A good strategy is to have as long and
"representative" sequence as possible as reference.

The reference sequence can be given in a separate file (see [Options](#options)
below), otherwise the first sequence in input will be used as reference.

If any sequence is reverse-complemented, the string `@revcomp@` is appended to
the fasta header.

If some sequences are too dissimilar to the reference sequence, they are
written to an (optional) file (see [Options](#options)).

## Usage

    $ reorient.pl infile.fasta
    $ reorient.pl --help

## Options

* `-r, --ref=<ref.fas>` Use (first) sequence in file ref.fas as reference. If
  ref.fas is not given, the first sequence in input will be used as reference.
* `-o, --out=<out.fas>` Write to file out.fas instead of STDOUT.
* `-f, --fail=<fail.fas>` Sequences with 'No hit' with ref sequence will be
  stored in fail.fas.
* `-v, --version` Show version number.
* `-V, --verbose` Be verbose about input/output. `--noverbose` will turn off
  extra printing.

## Dependencies

The script uses perl, perldoc and perl module StandAloneBlastPlus.pm. In
addition, ncbi-blast+ ("stand-alone blast") needs to be installed.

On a deb-based Linux system, they can be installed using:

    $ sudo apt install ncbi-blast+ bioperl bioperl-run

## Notes

Reverse-complemented sequences in the output can easily be identified using,
e.g., grep: `grep '@revcomp@' out.fas`.

The added string can be removed using GNU sed:
`sed -i 's/@revcomp@$//' out.fas`.
Or in one go:

    $ ./src/reorient.pl data/infile.fasta | \
        sed 's/@revcomp@$//' > out.fas

The success of identifying the sequence as being "plus" or "minus" depends on
the success of blastn finding a significant hit. This may sometimes not be
possible.  One may wish, however, to pay extra attention to those sequences as
they probably will be difficult to align (if this is the aim). To run the
script while saving too dissimilar sequences to a separate file, use the `-f`
option. For example:

    $ ./src/reorient.pl -f fail.fas data/infile.fasta

The functionality provided in this script is already implemented in, e.g.,
[`orient` command in
USEARCH](https://drive5.com/usearch/manual/cmd_orient.html), and option
[`--adjustdirection` in
MAFFT](https://mafft.cbrc.jp/alignment/software/adjustdirection.html).  Both
these examples uses another criterion for deciding the sequence direction.

