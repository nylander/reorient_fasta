# re-orient fasta sequences

- Last modified: tis jan 16, 2024  04:10
- Sign: nylander

## Description

Try to correct the direction of DNA sequences in fasta (by
reverse complement) in relation to a reference sequence.

More specifically, try to "re-orient" (reverse+complement) sequences in
infile.fasta based on their blast hits (bl2seq, blastn).  All sequences in
input are blasted against the same reference sequence. A good strategy is to
have as long and "representative" sequence as possible as reference.

The reference sequence can be given in a separate file (see options below),
otherwise the first sequence in input will be used as reference.

If any sequence is reverse-complemented, the string `@revcomp@` is appended to
the fasta header.  These sequences can then easily be identified by, e.g., using
grep: `grep '@revcomp@' file.fas`, or be removed by, e.g., using GNU sed:
`sed -i 's/@revcomp@$//' file.fas`.

## Usage

    $ reorient.pl infile.fasta
    $ reorient.pl --help

## Options

* `-r, --ref=<ref.fas>` Use (first) sequence in file ref.fas as reference. If
  ref.fas is not given, the first sequence in input will be used as reference.
* `-o, --out=<out.fas>` Write to out.fas instead of STDOUT.
* `-f, --fail=<fail.fas>` Sequences with 'No hit' with ref sequence will be
  stored in fail.fas.
* `-v, --version` Show version number
* `-V, --verbose` Be verbose about input/output. `--noverbose` will turn off
  extra printing.

## Dependencies

The script uses perl with perldoc and perl modules BioPerl and bl2seq
(which requires stand-alone blast).

## Notes

This functionality is already implemented in, e.g., [`orient` command in
USEARCH](https://drive5.com/usearch/manual/cmd_orient.html), and option
[`--adjustdirection` in MAFFT](https://mafft.cbrc.jp/alignment/software/).
