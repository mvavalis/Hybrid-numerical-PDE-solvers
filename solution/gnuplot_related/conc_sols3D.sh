#!/bin/bash
conc_sol=sol_3d.gnuplot
conc_sol_small=sol_3d_small.gnuplot
rm $conc_sol
rm $conc_sol_small


for file in *
do
  if [ -f "$file" ]; then
    if echo "$file" | grep -q 'solution-3d'.*'.gnuplot$'; then
      grep -v '#' $file >> $conc_sol
    fi
  fi
done

grep '#\$' $file | uniq >> $conc_sol


if [ -z "$1" ]; then
  sol_size=`wc -l $conc_sol | cut -f1 -d\ `
  better_size=20000
  kill_frac=$(($(($sol_size))/$(($better_size))))
  batch_kill=4  #nof lines
  
  cnt=0
  loop_bound=$(($kill_frac*$batch_kill))
  while read line
  do
    if [ $cnt -lt $batch_kill ]; then
      echo $line >> $conc_sol_small
    fi
    
    cnt=$(($cnt + 1))
    if [ $cnt -eq $loop_bound ]; then
      cnt=0
    fi
  done < $conc_sol
  
  grep '#\$' $file | uniq >> $conc_sol_small
fi
