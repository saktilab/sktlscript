system=$1
/home/ala/opt/amber20/bin/antechamber -i ${system}.pdb -fi pdb -o ${system}.mol2 -fo mol2 -c bcc -s 2 -nc 0
/home/ala/opt/amber20/bin/parmchk2 -i ${system}.mol2 -f mol2 -o ${system}.frcmod
cp ${system}.mol2 output.mol2
