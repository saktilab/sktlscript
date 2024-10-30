#!/usr/bin/env python3

import sys
import numpy as np

traj = open("traject.xyz", "w")
vel = open("velocity.xyz", "w")
 

with open(sys.argv[1], "r") as f:
    step = -1
    for line in f:
        coords = []
        veloc = []
        nat = int(line)
        info = next(f)
        for i in range(0, nat):
            line = next(f)
            arr = line.split()
            coords.append(arr[0:4])
            veloc.append(arr[5:8])     
        step +=1
        
        traj.write(str(nat)) 
        traj.write('\n')
        traj.write(info)
        for coord in coords:
            traj.write(' '.join(map(str, coord)))
            traj.write('\n')
        vel.write(str(nat)) 
        vel.write('\n')
        vel.write(info)
        for velocity in veloc:
            if (sys.argv[2] == "fs"):
                velocity = [float(x)*0.001 for x in velocity]
                vel.write(' '.join(map(str, velocity)))
            else:
                vel.write(' '.join(map(str, velocity)))
            vel.write('\n')
            