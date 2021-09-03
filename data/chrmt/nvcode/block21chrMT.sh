for individual in NA21141 NA21142 NA21143 NA21144 
do
 nohup bcftools consensus -s $individual -i 'VT!="SV"' -f /mnt/tmpfs/human-chrMT.fasta /mnt/tmpfs/chrMTnormalization.vcf.gz > $individual.fa 2> $individual.log
done
