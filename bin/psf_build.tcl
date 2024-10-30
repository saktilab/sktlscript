package require psfgen
topology water/water.par
topology oh/oh.par
topology 3qcdp/3qcdp.par
guesscoord
coordpdb meamix.pdb A
writepdb meamix.pdb
writepsf meamix.psf
