#!/bin/bash


filename=device_dedispersion_kernel.cu
rm $filename

cat device_dedispersion_kernel_header_template.cu > $filename

for range in {0,1,2,3,4,5,6,7}
do
	cat device_dedispersion_kernel_template.cu | sed "s/{RANGE}/_$range/g" >> $filename
done

cat device_dedispersion_kernel_footer_template.cu >> $filename

