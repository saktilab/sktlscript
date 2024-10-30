#!/bin/bash
for i in `ls bandstr/*_* -d`
do 
    cd $i
    rm INCAR KPOINTS
    cp ../../run_vasp.csh .
    #sampling
    cp INCAR_sampling INCAR
    cp KPOINTS_sampling KPOINTS

    #echo "Running bandstr calculation (charge): $i" >> $rootdir/auto.log
    #mpirun -np $NSLOTS /home/users/ahq/software/vasp_5.4.4/bin/vasp_std > vasp_sampling.out

    #mv OUTCAR OUTCAR.sampling
    #mv vasprun.xml vasprun.xml.sampling
    #mv PROCAR PROCAR.sampling
    #mv DOSCAR DOSCAR.sampling

    #cp INCAR_bs INCAR
    #cp KPOINTS_bs KPOINTS

    #echo "Running bandstr calculation (kpath): $i" >> $rootdir/auto.log
    #mpirun -np $NSLOTS /home/users/ahq/software/vasp_5.4.4/bin/vasp_std > vasp_bs.out
  
    cd ../../
done 

