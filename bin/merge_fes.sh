for i in 1 2;
do
 mkdir fes 
 cp ../run${i}/fes.dat fes/fes_${i}
 lines=$(sed -n '/500 GA/,/501 GA/p' fes.dat | wc -l)
 nlines=$(echo `expr $lines - 1`)
 var1=$(echo `expr $nlines + 2`)
 var2=$(echo `expr $nlines + 1`)
 ls -v fes/* | xargs -n 1 tail -n +${var1} > x
done
 head -n ${var2} fes/fes_1 > header
 cat header x > fes.dat
 rm x
 rm header
rm -rf fes
