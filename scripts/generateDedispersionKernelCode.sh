#!/bin/bash

maxRange=${1-7}

kernelFilename=../lib/device_dedispersion_kernel.cu
kernelTemplatePrefix=../lib/device_dedispersion_kernel
rm $kernelFilename

cat "$kernelTemplatePrefix"_header_template.cu > $kernelFilename

for range in `seq 0 $maxRange`
do
	cat "$kernelTemplatePrefix"_template.cu | sed "s/__RANGE__/_$range/g" >> $kernelFilename
done

cat "$kernelTemplatePrefix"_footer_template.cu >> $kernelFilename


kernelFunctionsFilename=../lib/headers/kernel_functions.h
kernelFunctionsTemplatePrefix=../lib/headers/kernel_functions
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

cat "$kernelFunctionsTemplatePrefix"_template.h | sed "s/{UNROLLS_RANGES}/$unrolls/g" | sed "s/{SNUMREG_RANGES}/$snumreg/g" | sed "s/{SDIVINT_RANGES}/$sdivint/g" | sed "s/{SDIVINDM_RANGES}/$sdivindm/g" | sed "s/{KERNEL_FUNCTION_RANGES}/$kernelFunctions/g" >> $kernelFunctionsFilename
