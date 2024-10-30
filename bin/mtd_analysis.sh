cp ../simu01/coord.inp .
cp ../simu01/3H.out .
cp ../simu01/gro.py .
cp ../simu01/list_nac.py .
./3H.out < coord.inp
./list_nac.py -e O -l dftb.inp -parm 5 traject mulliken 
./gro.py charge_list.dat -l dftb.inp
