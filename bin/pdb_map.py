#!/usr/bin/env python3
import sys

chaintype = []
atom_id = []
atom_type = []
resname = []
chain_id = []
resid = []
x = []
y = []
z = []
beta1 = []
beta2 = []
atom_name = []

with open(sys.argv[1], "r") as f:
    header = next(f)
    for line in f:
        if "END" in line:
            break
        arr = line.split()
        chaintype.append(arr[0])
        atom_id.append(int(arr[1])) 
        atom_type.append(arr[2])
        resname.append(arr[3])
        chain_id.append(arr[4])
        resid.append(int(arr[5]))
        x.append(float(arr[6]))
        y.append(float(arr[7]))
        z.append(float(arr[8]))
        beta1.append(float(arr[9]))
        beta2.append(float(arr[10]))
        atom_name.append(arr[11])

#START MAPPING
index = [i for i, word in enumerate(resname) if word.endswith("GLU")]
for i in index:
    if atom_type[i] == "HG2":
        atom_type[i] = "HG3" 
    if atom_type[i] == "HG1":
        atom_type[i] = "HG2"
    if atom_type[i] == "HB2":
        atom_type[i] = "HB3"
    if atom_type[i] == "HB1":
        atom_type[i] = "HB2"

index = [i for i, word in enumerate(resname) if word.endswith("ASP")]
for i in index:
    if atom_type[i] == "HB2":
        atom_type[i] = "HB3"
    if atom_type[i] == "HB1":
        atom_type[i] = "HB2"

index = [i for i, word in enumerate(resname) if word.endswith("LEU")]
for i in index:
    if atom_type[i] == "HB2":
        atom_type[i] = "HB3"
    if atom_type[i] == "HB1":
        atom_type[i] = "HB2"
    if atom_type[i] == "1HD1":
        atom_type[i] = "HD11"
    if atom_type[i] == "2HD1":
        atom_type[i] = "HD12"
    if atom_type[i] == "3HD1":
        atom_type[i] = "HD13"
    if atom_type[i] == "1HD2":
        atom_type[i] = "HD21"
    if atom_type[i] == "2HD2":
        atom_type[i] = "HD22"
    if atom_type[i] == "3HD2":
        atom_type[i] = "HD23"

index = [i for i, word in enumerate(resname) if word.endswith("LEU")]
for i in index:
    if atom_type[i] == "1HG1":
        atom_type[i] = "HG11"
    if atom_type[i] == "2HG1":
        atom_type[i] = "HG12"
    if atom_type[i] == "3HG1":
        atom_type[i] = "HG13"
    if atom_type[i] == "1HG2":
        atom_type[i] = "HG21"
    if atom_type[i] == "2HG2":
        atom_type[i] = "HG22"
    if atom_type[i] == "3HG2":
        atom_type[i] = "HG23"
    
index = [i for i, word in enumerate(resname) if word.endswith("THR")]
for i in index:
    if atom_type[i] == "1HG2":
        atom_type[i] = "HG21"
    if atom_type[i] == "2HG2":
        atom_type[i] = "HG22"
    if atom_type[i] == "3HG2":
        atom_type[i] = "HG23"

index = [i for i, word in enumerate(resname) if word.endswith("ARG")]
for i in index:
    if atom_type[i] == "HB2":
        atom_type[i] = "HB3"
    if atom_type[i] == "HB1":
        atom_type[i] = "HB2"
    if atom_type[i] == "HG2":
        atom_type[i] = "HG3"
    if atom_type[i] == "HG1":
        atom_type[i] = "HG2"
    if atom_type[i] == "HD2":
        atom_type[i] = "HD3"
    if atom_type[i] == "HD1":
        atom_type[i] = "HD2"
    if atom_type[i] == "2HH1":
        atom_type[i] = "HH12"
    if atom_type[i] == "1HH2":
        atom_type[i] = "HH21"
    if atom_type[i] == "2HH2":
        atom_type[i] = "HH22"

index = [i for i, word in enumerate(resname) if word.endswith("GLY")]
for i in index:
    if atom_type[i] == "HA2":
        atom_type[i] = "HA3"
    if atom_type[i] == "HA1":
        atom_type[i] = "HA2"

index = [i for i, word in enumerate(resname) if word.endswith("ILE")]
for i in index:
    if atom_type[i] == "CD":
        atom_type[i] = "CD1"
    if atom_type[i] == "2HH1":
        atom_type[i] = "HH12"
    if atom_type[i] == "HD1":
        atom_type[i] = "HD11"
    if atom_type[i] == "HD2":
        atom_type[i] = "HD12"
    if atom_type[i] == "HD3":
        atom_type[i] = "HD13"
    if atom_type[i] == "2HH1":
        atom_type[i] = "HH12"
    if atom_type[i] == "1HG2":
        atom_type[i] = "HG21"
    if atom_type[i] == "2HG2":
        atom_type[i] = "HG22"
    if atom_type[i] == "3HG2":
        atom_type[i] = "HG23"
    if atom_type[i] == "1HG1":
        atom_type[i] = "HG12"
    if atom_type[i] == "2HG1":
        atom_type[i] = "HG13"    

index = [i for i, word in enumerate(resname) if word.endswith("LYS")]
for i in index:
    if atom_type[i] == "HB2":
        atom_type[i] = "HB3"
    if atom_type[i] == "HB1":
        atom_type[i] = "HB2"    
    if atom_type[i] == "HG2":
        atom_type[i] = "HG3"
    if atom_type[i] == "HG1":
        atom_type[i] = "HG2"    
    if atom_type[i] == "HD2":
        atom_type[i] = "HD3"
    if atom_type[i] == "HD1":
        atom_type[i] = "HD2"    
    if atom_type[i] == "HE2":
        atom_type[i] = "HE3"
    if atom_type[i] == "HE1":
        atom_type[i] = "HE2"    
    
index = [i for i, word in enumerate(resname) if word.endswith("PRO")]
for i in index:
    if atom_type[i] == "HD2":
        atom_type[i] = "HD3"
    if atom_type[i] == "HD1":
        atom_type[i] = "HD2"    
    if atom_type[i] == "HG2":
        atom_type[i] = "HG3"
    if atom_type[i] == "HG1":
        atom_type[i] = "HG2"    
    if atom_type[i] == "HB2":
        atom_type[i] = "HB3"
    if atom_type[i] == "HB1":
        atom_type[i] = "HB2"    
    
index = [i for i, word in enumerate(resname) if word.endswith("SER")]
for i in index:
    if atom_type[i] == "HB2":
        atom_type[i] = "HB3"
    if atom_type[i] == "HB1":
        atom_type[i] = "HB2"    

index = [i for i, word in enumerate(resname) if word.endswith("PHE")]
for i in index:
    if atom_type[i] == "HB2":
        atom_type[i] = "HB3"
    if atom_type[i] == "HB1":
        atom_type[i] = "HB2"    

index = [i for i, word in enumerate(resname) if word.endswith("MET")]
for i in index:
    if atom_type[i] == "HB2":
        atom_type[i] = "HB3"
    if atom_type[i] == "HB1":
        atom_type[i] = "HB2"    
    if atom_type[i] == "HG2":
        atom_type[i] = "HG3"
    if atom_type[i] == "HG1":
        atom_type[i] = "HG2"    
    
index = [i for i, word in enumerate(resname) if word.endswith("GLN")]
for i in index:
    if atom_type[i] == "HB2":
        atom_type[i] = "HB3"
    if atom_type[i] == "HB1":
        atom_type[i] = "HB2"   
    if atom_type[i] == "HG2":
        atom_type[i] = "HG3"
    if atom_type[i] == "HG1":
        atom_type[i] = "HG2"   
    if atom_type[i] == "1HE2":
        atom_type[i] = "HE21"
    if atom_type[i] == "2HE2":
        atom_type[i] = "HE22"   

index = [i for i, word in enumerate(resname) if word.endswith("TRP")]
for i in index:
    if atom_type[i] == "HB2":
        atom_type[i] = "HB3"
    if atom_type[i] == "HB1":
        atom_type[i] = "HB2"  

index = [i for i, word in enumerate(resname) if word.endswith("CYS")]
for i in index:
    if atom_type[i] == "HB2":
        atom_type[i] = "HB3"
    if atom_type[i] == "HB1":
        atom_type[i] = "HB2"  

index = [i for i, word in enumerate(resname) if word.endswith("TYR")]
for i in index:
    if atom_type[i] == "HB2":
        atom_type[i] = "HB3"
    if atom_type[i] == "HB1":
        atom_type[i] = "HB2"  

index = [i for i, word in enumerate(resname) if word.endswith("ASN")]
for i in index:
    if atom_type[i] == "HB2":
        atom_type[i] = "HB3"
    if atom_type[i] == "HB1":
        atom_type[i] = "HB2"  
    if atom_type[i] == "1HD2":
        atom_type[i] = "HD21"
    if atom_type[i] == "2HD2":
        atom_type[i] = "HD22"  

index = [i for i, word in enumerate(resname) if word.endswith("HIS")]
for i in index:
    if atom_type[i] == "HB2":
        atom_type[i] = "HB3"
    if atom_type[i] == "HB1":
        atom_type[i] = "HB2"  
    
index = [i for i, word in enumerate(resname) if word.endswith("HIP")]
for i in index:
    if atom_type[i] == "HB2":
        atom_type[i] = "HB3"
    if atom_type[i] == "HB1":
        atom_type[i] = "HB2"  


for chaintype, atom_id, atom_type, resname, chain_id, resid,  x, y, z, beta1, beta2, atom_name in zip(chaintype, atom_id, atom_type, resname, chain_id, resid, x, y, z, beta1, beta2, atom_name):
    print('{}  {:5d} {:>4s} {} {} {:3d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}      {}'.format(chaintype, atom_id, atom_type, resname, chain_id, resid, x, y, z, beta1, beta2, atom_name))
print("END")    