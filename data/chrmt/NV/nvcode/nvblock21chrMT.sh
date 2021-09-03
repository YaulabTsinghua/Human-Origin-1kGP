for individual in NA21141 NA21142 NA21143 NA21144 
do
 nohup perl NVdna.pl $individual.fa > $individual.log 2>&1
done
