#!/bin/csh
#
set echo
#
cp spline.H /$HOME/include
#
g++ -c -g spline.C >& compiler.out
if ( $status != 0 ) then
  echo "Errors compiling spline.C."
  exit
endif
rm compiler.out
#
mv spline.o ~/libcpp/$ARCH/spline.o
#
echo "A new version of spline has been created."
