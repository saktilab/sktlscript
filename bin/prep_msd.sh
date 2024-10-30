for i in $@
do
  mkdir ${i} 
  (
  cd ${i} || continue
  rm dftb.traj dftb.out dftb.inp
  ln -s ~ahq/work/Li/qmd/solution/${i}/nve/dftb.traj dftb.traj
  ln -s ~ahq/work/Li/qmd/solution/${i}/nve/dftb.inp dftb.inp
  ln -s ~ahq/work/Li/qmd/solution/${i}/nve/dftb.out dftb.out
  python3.4 ~/bin/make_diffusion_input.py ${i} > diff.inp
  )
done
