file mkdir psf
exec rm -f psf/*pdb psf/*psf

#Base name of your pdb & psf files
set name ice
set segment_name TIP4

#Preparing your input pdb file
#In many cases, you may need to customize the preparation
exec cp $name.pdb psf/$name.pdb

#Embed the psfgen command in this script
#psfgen << ENDMOL
package require psfgen
resetpsf
#Insert your topology file
topology tip4p.top

#Build your pdb segment
segment $segment_name {
 pdb psf/$name.pdb
}


#Read the coordinate
coordpdb psf/$name.pdb $segment_name

#Guess missing coordinate
guesscoord

#write structure and coordinate files
writepsf psf/$name.psf
writepdb psf/$name.pdb

