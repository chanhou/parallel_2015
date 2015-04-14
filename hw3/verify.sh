#!/bin/bash
src="/home/course/hw03"
a='test1.dat test2.dat'
for i in $a; do
	# echo $i
	./mainAOS.exe "$src/$i" "$i-AOS.out"
	./mainSOA.exe "$src/$i" "$i-SOA.out"
	$src/mainREF.exe "$src/$i" "$i-REF.out"
	a=`diff $i-AOS.out $i-REF.out`
	b=`diff $i-SOA.out $i-REF.out`

	if [ -z "$a" ]; then
		echo "mainAOS.exe verified successfully with $i dataset"
	fi
	if [ -z "$b" ]; then
		echo "mainSOA.exe verified successfully with $i dataset"
	fi
done
