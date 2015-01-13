#!/usr/bin/perl 

use strict;

my $gp=q{
set terminal pdf 
set out "#file#.pdf"
set xlabel "# threads"
set ylabel "#y#"
show xlabel 
show ylabel
set yrange [0:]
set xrange [0:]
plot "__1.dat" using 1:2 with linespoints pt 7 , \
};

open(D1,">","__1.dat");
open(D2,">","__2.dat");
my $base;

while(<>){
  print D1 $_;

  if   (/^0\s+(\d+(.\d*)?)/){ 
      $base=$1 ; 
      print D2 "0 1\n"}
  elsif(/^(\d+)\s+(\d+(.\d*)?)/){ 
      print D2 "$1 ", ($base/$2),"\n"}
}
close D1;
close D2;

open(GP,"|-","gnuplot" );
print GP $gp =~ s/#y#/time/r =~ s/#file#/time/r;;
close GP;

open(GP,"|-","gnuplot" );
$gp =~ s/1\./2./g;
print GP $gp =~ s/#y#/speedup/r =~ s/#file#/speed/r;
close GP;

