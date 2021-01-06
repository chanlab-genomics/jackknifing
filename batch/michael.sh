#!/bin/bash
#PBS -l select=1:ncpus=8:mem=32GB
#PBS -l walltime=0:05:00
#PBS -A NCMAS-d85
#PBS -k oe
#PBS -j oe
#PBS -N BLASTN_self

cd $PBS_O_WORKDIR
module load blast
for file in $(ls ./Sequence_Files/*.fa)
do
        file_basename=$(basename $file | cut -f 1 -d '.')
        blastn -query $file  -subject $file -outfmt '5' -out Self_Blasts/${file_basename}_self_blast.xml
        xsltproc --novalid ~/Programs/xslt-sandbox-master/xslt-sandbox-master/stylesheets/bio/ncbi/blast2dotplot.xsl Self_Blasts/${file_basename}_self_blast.xml > ${file_basename}.html
done