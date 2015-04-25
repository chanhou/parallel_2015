#!/bin/bash
rm result_AOS.out result_SOA.out

for i in {1..6}
do
   # echo "Welcome $i times"
   # ./generate $i 1000 3 10 "test-$i-dim.dat"
   ./mainAOS.exe "test-$i-dim.dat" "test-$i-dim-AOS.out" >> result_AOS.out
   ./mainSOA.exe "test-$i-dim.dat" "test-$i-dim-SOA.out" >> result_SOA.out
done