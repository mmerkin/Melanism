#!/bin/bash

REF=Pp_scaffolds.fa
THREADS=32
datapath=/pub64/mattm/SV_calling/Pp_insurveyor/Pp_samples
outdir=/pub64/mattm/SV_calling/Pp_insurveyor/merge_vcfs

for filepath in "$datapath"/*; do
filename=${filepath##*/}
cp $filepath/out.pass.vcf.gz $outdir/${filename}_inserts.vcf.gz
done

echo "completed"

