#!/bin/bash
src="."
a='test1.dat test2.dat'
# make
for i in $a; do
	echo $i
	OMP_NUM_THREADS=8 ./mainAOSp.exe "$src/$i" "$i-AOS.out"
	OMP_NUM_THREADS=8 ./mainSOAp.exe "$src/$i" "$i-SOA.out"
	# $src/mainREF.exe "$src/$i" "$i-REF.out"
	a=`diff $i-AOS.out $i-REF.out`
	b=`diff $i-SOA.out $i-REF.out`

	if [ -z "$a" ]; then
		echo "mainAOS.exe verified successfully with $i dataset"
	fi
	if [ -z "$b" ]; then
		echo "mainSOA.exe verified successfully with $i dataset"
	fi
done
