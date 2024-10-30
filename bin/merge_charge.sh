mkdir charge 
for i in $@
do 
cp ../run${i}/mulliken charge/mull_${i}
done
natom=$(head -n 1 charge/mull_1)
var1=$(echo `expr $natom + 3`)
var2=$(echo `expr $natom + 2`)
ls -v charge/* | xargs -n 1 tail -n +${var1} > x
head -n ${var2} charge/mull_1 > head
cat head x > mulliken
rm x
rm head
rm -rf charge
