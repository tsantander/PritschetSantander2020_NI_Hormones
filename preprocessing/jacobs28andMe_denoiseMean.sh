#!/bin/bash

# Denoise and bias-regularise mean EPIs.

cd /Volumes/LTD/28andMe/allData

for d in ./*; do

	thisDir="${d%*/}"
	
	inFile=($(find "${thisDir}/data.functional.mean" -maxdepth 1 -type f -name "meanuf2018*.nii"))
	
	DenoiseImage -d 3 \
		-i "${inFile}" \
		-o "${thisDir}/data.functional.mean/denoisedMean.nii.gz" \
		-v 1
	
	N3BiasFieldCorrection 3 \
		"${thisDir}/data.functional.mean/denoisedMean.nii.gz" \
		"${thisDir}/data.functional.mean/meanN3.nii.gz"
		
done