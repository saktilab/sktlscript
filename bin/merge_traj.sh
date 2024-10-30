mkdir trajectory
for i in $@
do 
cp ../run${i}/dftb.traj trajectory/traj_${i}
done
natom=$(head -n 1 trajectory/traj_1)
var1=$(echo `expr $natom + 3`)
var2=$(echo `expr $natom + 2`)
ls -v trajectory/* | xargs -n 1 tail -n +${var1} > x
head -n ${var2} trajectory/traj_1 > head
cat head x > dftb.traj
rm x
rm head
rm -rf trajectory
