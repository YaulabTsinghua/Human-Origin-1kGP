for individual in NA21141 NA21142 NA21143 NA21144 
do
 nohup perl NVdna_20190908.pl $individual.fa 4 > $individual.log 2>&1
done
