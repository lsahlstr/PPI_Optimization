#!/usr/bin/env perl

use strict;
use warnings;

my $res1 = $ARGV[0];    # last residue in chain A
my $res2 = $ARGV[1];    # last residue in chain B
my $nframes = $ARGV[2]; # number of snapshots

my $string = "$res1 $res2";

open(IN,"pairdist.out");

# Split pairdist.out into a separate file for each frame
my $c = 1; # counter for number of snapshots
`mkdir -p pairdist`;
open(OUT,">pairdist/pairdist.$c");

while (my $l = <IN>) {
    if ($l !~ /^$string/) {
        print OUT $l;
    } elsif ($l =~ /^$string/) {
        print OUT $l;
        close OUT;
        ++$c;
        if ($c <= $nframes) {
            open(OUT,">pairdist/pairdist.$c");
        }
    }
}
