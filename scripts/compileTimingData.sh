#!/bin/bash

rm timing.csv

echo "range, Speedup with Range Specific Params, errorHigh, errorLow" >> timing.csv
for range in {0,1,2,3,4,5,6,7}
do
	cat range${range}Stats.txt | tail -1 | awk -v var="$range" -F " " '{print var", "$5", "$7-$5", "$5-$9}' >> timing.csv
done

paramString=$(cat stats.txt | tail -1 | awk -F " " '{print $3}')

echo >> timing.csv
echo "range, Speedup with avg Params, errorHigh, errorLow" >> timing.csv
for range in {0,1,2,3,4,5,6,7}
do
	cat range${range}Stats.txt | grep "$paramString" | awk -v var="$range" -F " " '{print var", "$4", "$6-$4", "$4-$8}' >> timing.csv
done

echo >> timing.csv
echo "range, UNROLLS, SNUMREG, SDIVINT, SDIVINDM" >> timing.csv
echo $paramString | sed 's/[^0-9^_]//g' | awk -F "_" '{print "avg, "$1", "$2", "$3", "$4}' >> timing.csv

for range in {0,1,2,3,4,5,6,7}
do
	cat range${range}Stats.txt | tail -1 | awk -F " " '{print $3}' | sed 's/[^0-9^_]//g' | awk -v var="$range" -F "_" '{print var", "$1", "$2", "$3", "$4}' >> timing.csv
done







