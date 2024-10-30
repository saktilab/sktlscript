moveatoms=$(cat move_atoms.dat)
frozenatoms=$(cat frozen_atoms.dat)
for i in $frozenatoms;
do gawk -vfroze="$i" 'NR==froze {print $1,$2,$3,$4" *"}' $1 > line_${i};
done

for i in $moveatoms;
do gawk -vmove="$i" 'NR==move {print $1,$2,$3,$4}' $1 > line_${i};
done

cat line_* > geom.xyz; 
rm line_*
