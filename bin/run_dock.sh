mkdir 2_gen_sphere 3_gen_grid 4_docking
ln -s -r /home/aditya/opt/dock6/parameters/

#Entering 2_gen_shphere to generate the spheres
cd 2_gen_sphere
echo "Generating spheres"
dms ../1_prep/cutted_noH.pdb -n -w 1.4 -v -o pro.ms
cp /home/aditya/Joseph/dock/gen_sphere/INSPH .
~/opt/dock6/bin/sphgen
~/opt/dock6/bin/sphere_selector pro.sph ../1_prep/ligand.mol2 10.0

#Creating box
cd ../3_gen_grid
echo "Generating grid"
cp /home/aditya/Joseph/dock/gen_grid/box.in .
~/opt/dock6/bin/showbox < box.in
#Creating grid
cp /home/aditya/Joseph/dock/gen_grid/grid.in .
~/opt/dock6/bin/grid -i grid.in

#Creating NChem grid Zou's Generalized Born
#mkdir nchemgrid_GB
#cd nchemgrid_GB
#cp /home/aditya/opt/dock6/tutorials/solvent_scoring_demo/3_grid/nchemgrid_GB/cavity.pdb .
#cp /home/aditya/opt/dock6/tutorials/solvent_scoring_demo/3_grid/nchemgrid_GB/INCHEM .
#echo "Generating grid for Zou's generalized Born model"
#~/opt/dock6/bin/nchemgrid_GB 
#cd ../
#Creating NChem grid Zou's Surface Area
#mkdir nchemgrid_SA
#cd nchemgrid_SA
#cp /home/aditya/opt/dock6/tutorials/solvent_scoring_demo/3_grid/nchemgrid_SA/INCHEM .
#echo "Generating grid for Zou's surface area model"
#~/opt/dock6/bin/nchemgrid_SA


#Performing docking
echo "Performing Docking"
cd ../4_docking
cp /home/aditya/Joseph/dock/docking/flex.in .
echo "Performing primary score flexible docking"
cp ../3_gen_grid/grid.* .
~/opt/dock6/bin/dock6 -i flex.in -o flex.out 
#echo "Performing GBSA Zou"
#cp /home/aditya/opt/dock6/tutorials/solvent_scoring_demo/4_dock/gbsa_zou.in .
#~/opt/dock6/bin/dock6 -i gbsa_zou.in -o gbsa_zou.out 
#echo "Performing GBSA Hawkins"
#cp /home/aditya/opt/dock6/tutorials/solvent_scoring_demo/4_dock/gbsa_hawkins.in .
#~/opt/dock6/bin/dock6 -i gbsa_hawkins.in -o gbsa_hawkins.out 

#DONE


