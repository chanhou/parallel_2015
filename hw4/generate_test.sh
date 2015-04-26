#!/bin/bash
# rm result_AOS.out result_SOA.out

for i in {1..6}
do
   # echo "Welcome $i times"
   ./generate $i 1000 3 10 "test-$i-dim.dat"
   ./mainAOS.exe "test-$i-dim.dat" "test-$i-dim-AOS.out" >> result_AOS_dim.out
   ./mainSOA.exe "test-$i-dim.dat" "test-$i-dim-SOA.out" >> result_SOA_dim.out
   ./mainAOSp.exe "test-$i-dim.dat" "test-$i-dim-AOSp.out" >> result_AOSp_dim.out
   ./mainSOAp.exe "test-$i-dim.dat" "test-$i-dim-SOAp.out" >> result_SOAp_dim.out
done

for i in {1..6}
do
   # echo "Welcome $i times"
   let "z=10**$i"
   echo $z
   ./generate 2 $z 3 10 "test-$i-point.dat"
   ./mainAOS.exe "test-$i-point.dat" "test-$z-point-AOS.out" >> result_AOS_point.out
   ./mainAOSp.exe "test-$i-point.dat" "test-$z-point-AOSp.out" >> result_AOSp_point.out
   ./mainSOA.exe "test-$i-point.dat" "test-$z-point-SOA.out" >> result_SOA_point.out
   ./mainSOAp.exe "test-$i-point.dat" "test-$z-point-SOAp.out" >> result_SOAp_point.out
done

for i in {1..6}
do
   # echo "Welcome $i times"
   let "z=10**$i"
   echo $z
   ./generate 2 1000 3 $z "test-$i-grid.dat"
   ./mainAOS.exe "test-$i-grid.dat" "test-$z-grid-AOS.out" >> result_AOS_grid.out
   ./mainAOSp.exe "test-$i-grid.dat" "test-$z-grid-AOSp.out" >> result_AOSp_grid.out
   ./mainSOA.exe "test-$i-grid.dat" "test-$z-grid-SOA.out" >> result_SOA_grid.out
   ./mainSOAp.exe "test-$i-grid.dat" "test-$z-grid-SOAp.out" >> result_SOAp_grid.out
done