#!/bin/bash
~/miniconda3/bin/antechamber -i $1.pdb -fi pdb -o $1.mol2 -fo mol2 -c bcc -s 2 -nc $2
~/miniconda3/bin/parmchk2 -i $1.mol2 -f mol2 -o $1.frcmod
touch tleap.in
echo "source leaprc.gaff" >> tleap.in
echo "$1 = loadmol2 $1.mol2" >> tleap.in
echo "check $1" >> tleap.in
echo "loadamberparams $1.frcmod" >> tleap.in
echo "saveoff $1 $1.lib" >> tleap.in
echo "saveamberparm $1 $1.prmtop $1.inpcrd" >> tleap.in
echo "savepdb $1 $1_gaff.pdb" >> tleap.in
echo "quit" >> tleap.in
~/miniconda3/bin/tleap -f tleap.in
# acpype -p $1.prmtop -x $1.inpcrd
# resname=$(head -n 2 $1.mol2 | tail -n 1)
# gsed "s/$resname/$1/g" *.top > $1.top
