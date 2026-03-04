#!/usr/bin/env perl
use strict;
use warnings;

my $state = 0;
my $counter = 0;

while(<>) {
    if ($state == 1) { # $counter is fixed number of lines left to print
        print $_;
        if (--$counter == 0) {
            $state = 0;
        }
        next;
    }
    elsif ($state == 2) { # variable-length GPNORM block; expect line $counter or AVE if zero
        if ( $counter == 0 ? m{^\s*AVE\s} : m{^\s*$counter\s} ) {
            print $_;
            ++$counter;
            next;
        }
        else {
            $state = 0; $counter = 0; # block finished
        }
    }
    elsif ($state == 3) { # variable-length FULL-POS SPNORM block
        if ( m{^\s*(?:[A-Z][0-9]+\s+)?[0-9]+\s+:\s+} ) {
            print $_;
            ++$counter;
            next;
        }
        else {
            $state = 0; $counter = 0; # block finished
        }
    }
    elsif ($state == 4) { # variable-length FULL-POS GPNORM block; expect header if $counter zero
        if ( $counter == 0 ? m{^\s+AVERAGE\s+MINIMUM\s+MAXIMUM\b} : m{^\s*[0-9]+\s+:\s+} ) {
            print $_;
            ++$counter;
            next;
        }
        else {
            $state = 0; $counter = 0; # block finished
        }
    }
    
    # base state
    if ( m{^\s*SPECTRAL NORMS\s} ) {
        $state = 1; $counter = 2; print $_; # print with next two lines
    }
    elsif ( m{^\s*GPNORM\b} ) {
        $state = 2; $counter = 0; print $_; # print variable-length GPNORM block
    }
    elsif ( m{^\s*FULL-POS SPNORMS\b} ) {
        $state = 3; $counter = 0; print $_; # print variable-length FULL-POS SPNORM block
    }
    elsif ( m{^\s*FULL-POS GPNORMS\b} ) {
        $state = 4; $counter = 0; print $_; # print variable-length FULL-POS GPNORM block
    }
    elsif ( m{\sSP NORMS\s} ) {
        print $_; # print just this line
    }
    elsif ( m{^\s*MASSDIA[0-9]*\s} ) {
        print $_; # print just this line
    }
}
