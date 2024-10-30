mkdir 2_gen_sphere 3_gen_grid 4_docking

#Entering 2_gen_shphere to generate the spheres
cd 2_gen_sphere
cp ../1_prep/ligand.mol2 .
cp ../1_prep/cutted_noH.pdb .
dms cutted_noH.pdb -n -w 1.4 -v -o pro.ms
cp /home/aditya/Joseph/dock/gen_sphere/INSPH .
~/opt/dock6/bin/sphgen
~/opt/dock6/bin/sphere_selector pro.sph ligand.mol2 10.0

#Creating box
cd ../3_gen_grid
cp ../2_gen_sphere/selected_spheres.sph .
cp ../1_prep/ligand.mol2 .
cp /home/aditya/Joseph/dock/gen_grid/box.in .
~/opt/dock6/bin/showbox < box.in
#Creating grid
cp /home/aditya/Joseph/dock/gen_grid/grid.in .
~/opt/dock6/bin/grid -i grid.in

#Performing docking
cd ../4_docking
cp ../1_prep/ligand.mol2 .
cp ../2_gen_sphere/selected_spheres.sph .
cp ../3_gen_grid/grid.nrg .
cp ../3_gen_grid/grid.bmp .
cp /home/aditya/Joseph/dock/docking/flex.in .
~/opt/dock6/bin/dock6 -i flex.in -o flex.out 

#DONE


