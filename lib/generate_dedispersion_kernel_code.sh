#!/bin/bash

maxRange=${1-7}

kernelFilename=device_dedispersion_kernel.cu
rm $kernelFilename

cat device_dedispersion_kernel_header_template.cu > $kernelFilename

for range in `seq 0 $maxRange`
do
	cat device_dedispersion_kernel_template.cu | sed "s/{RANGE}/_$range/g" >> $kernelFilename
done

cat device_dedispersion_kernel_footer_template.cu >> $kernelFilename


kernelFunctionsFilename=AstroAccelerate/kernel_functions.h
rm $kernelFunctionsFilename

for range in `seq 0 $maxRange`
do
	unrolls+=UNROLLS_$range
	snumreg+=SNUMREG_$range
	sdivint+=SDIVINT_$range
	sdivindm+=SDIVINDM_$range
	kernelFunctions+=shared_dedisperse_kernel_range_$range
	if [ $range -lt $maxRange ]
	then
		unrolls+=,
		snumreg+=,
		sdivint+=,
		sdivindm+=,
		kernelFunctions+=,
	fi
done

cat AstroAccelerate/kernel_functions_template.h | sed "s/{UNROLLS_RANGES}/$unrolls/g" | sed "s/{SNUMREG_RANGES}/$snumreg/g" | sed "s/{SDIVINT_RANGES}/$sdivint/g" | sed "s/{SDIVINDM_RANGES}/$sdivindm/g" | sed "s/{KERNEL_FUNCTION_RANGES}/$kernelFunctions/g" >> $kernelFunctionsFilename
