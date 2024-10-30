for i in `seq 1 144`
do 
 sed -n '${i}p' $1 > $i.smi
done
