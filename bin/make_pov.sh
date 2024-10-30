~/bin/xyz2poscar.rb opt.xyz > POSCAR;
~/bin/poscar2pov.py -i POSCAR --rotate '1x,-30y,1z';
aflow --dist=3.0 < POSCAR 
