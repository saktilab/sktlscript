catdcd -o dftb.traj -otype xyz -stype pdb -s ../../prep/system_init.pdb -dcd system_out.dcd 
~/bin/diff_com < diff.inp
~/bin/md_msd_fit.py msd.txt > diff.txt
