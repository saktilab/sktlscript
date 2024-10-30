for i in $@
do
  (
 cd ${i}/output_nvt || continue
~solccp/bin/catdcd -otype xyz -o dftb.traj -stype pdb -s ../prep/system_init.pdb -stride 10 system_nvt.dcd
~/.pyenv/versions/3.5.1/bin/python ~/bin/make_diffusion_input.py ${i} > diff.inp
~/bin/diff_com < diff.inp 
~/.pyenv/versions/3.5.1/bin/python ~/bin/md_msd_fit.py -s 30 msd.txt > diff.txt
 )
done
