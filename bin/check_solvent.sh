for i in $@
do 
grep "number" ~solccp/work/data/Li-bat/CMD_GAFF/pure_solvent/${i}/prep/packmol.inp | awk '{print $2}'
done
