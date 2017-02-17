#!/bin/bash

rm ../input_files/tests/*
rm -rf ../output/57*
rm -rf ../output/tests/*

mkdir ../input_files/tests

for i in $(ls ~/filterbank/unit_tests/*)
do
	cat ../input_files/header > input.txt
	echo "file $i" >> input.txt
	j=$(echo $i | awk -F"/" '{print $6}')
	mv input.txt ../input_files/tests/$j.txt
done

for i in $(ls ../input_files/tests/*)
do
	echo $i
	./astro-accelerate.sh $i
done

for i in $(ls ../output/5*/gl*frb*)
do
	j=$(echo $i | awk -F"/" '{print $3}'| awk -F"." '{print $1}')
	echo $j
	echo "set terminal png" > $j.plt
	echo "set output \"$j.png\"" >> $j.plt
	echo "set title \"$j\"" >> $j.plt
	echo "splot \"./$i\" binary format=\"%float%float%float%float\" u 1:2:3 palette" >> $j.plt
	mv $j.plt ../output/tests/$j.plt
done

for i in $(ls ../output/tests/*.plt)
do
	gnuplot $i
done

convert -delay 150 -loop 0 *.png ani.gif

mv *.png ../output/tests/
mv *.gif ../output/tests/

google-chrome-stable ../output/tests/ani.gif &
