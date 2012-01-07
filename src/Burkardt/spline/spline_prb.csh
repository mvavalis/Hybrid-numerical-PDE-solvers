#!/bin/csh
#
set echo
#
g++ -c -g -I/$HOME/include spline_prb.C >& compiler.out
if ( $status != 0 ) then
  echo "Errors compiling spline_prb.C"
  exit
endif
rm compiler.out
#
g++ spline_prb.o /$HOME/libcpp/$ARCH/spline.o -lm
if ( $status != 0 ) then
  echo "Errors linking and loading spline_prb.o"
  exit
endif
#
rm spline_prb.o
#
mv a.out spline_prb
./spline_prb > spline_prb.out
if ( $status != 0 ) then
  echo "Errors running spline_prb"
  exit
endif
rm spline_prb
#
echo "The spline_prb program has been executed."
