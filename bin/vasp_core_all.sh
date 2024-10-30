#!/bin/bash
#PBS -l select=1:ncpus=36:mpiprocs=36:ompthreads=1:jobtype=core
#PBS -l walltime=120:00:00
#PBS -m be
#PBS -r n


SCRIPTS=/home/users/ahq/work/parameterization/scripts/

cd ${PBS_O_WORKDIR}

export PATH=$PATH:$SCRIPTS

NSLOTS=`cat $PBS_NODEFILE | wc -l`


#opt

rootdir=`pwd`

date > $rootdir/auto.log

echo "Optimizing the unit cell..." >> $rootdir/auto.log

cd opt

$SCRIPTS/make_vasp_kpoints.py

echo "Optimization 1" >> $rootdir/auto.log
mpirun -np $NSLOTS /home/users/ahq/software/vasp_5.4.4/bin/vasp_std > vasp.out


cp CONTCAR POSCAR
echo "Optimization 2" >> $rootdir/auto.log
mpirun -np $NSLOTS /home/users/ahq/software/vasp_5.4.4/bin/vasp_std > vasp.out

cp CONTCAR POSCAR
echo "Optimization 3" >> $rootdir/auto.log
mpirun -np $NSLOTS /home/users/ahq/software/vasp_5.4.4/bin/vasp_std > vasp.out

cp CONTCAR POSCAR
echo "Optimization 4" >> $rootdir/auto.log
mpirun -np $NSLOTS /home/users/ahq/software/vasp_5.4.4/bin/vasp_std > vasp.out

cp CONTCAR POSCAR
echo "Optimization 5" >> $rootdir/auto.log
mpirun -np $NSLOTS /home/users/ahq/software/vasp_5.4.4/bin/vasp_std > vasp.out


cd ..

echo "Setting up bandstr calculation" >> $rootdir/auto.log

$SCRIPTS/make_bandstr_input.py >> $rootdir/auto.log

for i in `ls bandstr/*_* -d`
do 
    cd $i
    rm INCAR KPOINTS

    #sampling
    cp INCAR_sampling INCAR
    cp KPOINTS_sampling KPOINTS

    echo "Running bandstr calculation (charge): $i" >> $rootdir/auto.log
    mpirun -np $NSLOTS /home/users/ahq/software/vasp_5.4.4/bin/vasp_std > vasp_sampling.out

    mv OUTCAR OUTCAR.sampling
    mv vasprun.xml vasprun.xml.sampling
    mv PROCAR PROCAR.sampling
    mv DOSCAR DOSCAR.sampling

    cp INCAR_bs INCAR
    cp KPOINTS_bs KPOINTS

    echo "Running bandstr calculation (kpath): $i" >> $rootdir/auto.log
    mpirun -np $NSLOTS /home/users/ahq/software/vasp_5.4.4/bin/vasp_std > vasp_bs.out
  
    cd ../../
done 


echo "Calculating the E-V curves" >> $rootdir/auto.log
$SCRIPTS/make_ev_input.py >> $rootdir/auto.log

cd e-v
for i in `ls *_* -d`
do
    cd $i
    echo "Running e-v static calculation: $i" >> $rootdir/auto.log
    mpirun -np $NSLOTS /home/users/ahq/software/vasp_5.4.4/bin/vasp_std > vasp.out
    cd ..;
done

$SCRIPTS/get_vasp_ev.py >> $rootdir/auto.log

cd results
for i in `ls * -d`
do
    cd $i
    $SCRIPTS/PLOT-birch_sol.py > output.log 
    cd ..
done
cd ..

cd ..


echo "Setting up elastic constant calculations" >> $rootdir/auto.log
$SCRIPTS/ElaStic_Setup_VASP >> $rootdir/auto.log


echo "Running elastic calculations" >> $rootdir/auto.log
cd elastic
for i in Dst0?/Dst??_??
do
    cd $i
    echo "Running static elastic calcuation: $i" >> $rootdir/auto.log
    mpirun -np $NSLOTS /home/users/ahq/software/vasp_5.4.4/bin/vasp_std > vasp.out
    cd ../../
done

echo "Analyzing and computing the results for elastic constants" >> $rootdir/auto.log

$SCRIPTS/ElaStic_Analyze_Energy   >> $rootdir/auto.log
$SCRIPTS/ElaStic_Result_Energy_2nd >> $rootdir/auto.log

cd ..


# surface energy











echo "The End" >> $rootdir/auto.log
date >> $rootdir/auto.log



