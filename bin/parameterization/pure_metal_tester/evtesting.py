#!/usr/bin/env python3

import numpy as np
import pymatgen as mg

def make_strain_structures(struct, dirn, mdr=0.05, NoP=11):
    scale = [] 
    if (scale==[]):
        ascale=1.
    else:
        ascale=scale[0]

    stretchstr = [] 
    if (stretchstr==[]):
        stretch=[1.,1.,1.]
    else:
        stretch=np.array(map(float,stretchstr[0].split()))

    bv = struct.lattice.matrix

    M_old= np.array(bv)
    D    = np.linalg.det(M_old)
    V0   = abs(stretch[0]*stretch[1]*stretch[2]*ascale**3*D)

    ptn = int((NoP-1)/2)
    structures = []

    cont= 0
    for s in range(-ptn, ptn+1):
        r = mdr*s/ptn
        if (s==0): r = 0.00001
        eps = r
        def_matrix={\
        'VOL'  :[[1.+eps      , 0.          , 0.             ],
                 [0.          , 1.+eps      , 0.             ],
                 [0.          , 0.          , 1+eps          ]],\

        'BOA'  :[[(1+eps)**-.5, 0.          , 0.             ],
                 [ 0.         , 1.+eps      , 0.             ],
                 [ 0.         , 0.          ,(1+eps)**-.5    ]],\

        'COA'  :[[(1+eps)**-.5, 0.          , 0.             ],
                 [ 0.         , (1+eps)**-.5, 0.             ],
                 [ 0.         , 0.          , 1.+eps         ]],\

        'ALPHA':[[1./(1-eps**2), 0.           , 0.           ],
                 [ 0.          , 1.           ,eps           ],
                 [ 0.          ,eps           , 1.           ]],\

        'BETA' :[[ 1.          , 0.           ,eps           ],
                 [ 0.          , 1./(1-eps**2), 0.           ],
                 [eps          , 0.           , 1.           ]],\

        'GAMMA':[[ 1.          ,eps           , 0.           ],
                 [eps          , 1.           , 0.           ],
                 [ 0.          , 0.           , 1./(1-eps**2)]]}
            
        M_eps = np.array(def_matrix[dirn])
        M_new = np.dot(M_old, M_eps)

        cont += 1
        # if (0  < cont and cont <  10): dirn_num = dirn.lower() + '_0'+ str(cont)
        # if (9  < cont and cont < 100): dirn_num = dirn.lower() + '_' + str(cont)

        newLattice = mg.Lattice(M_new)
                  

        # Writing the structure file ------------------------------------------------------------------
        newStr = mg.Structure(newLattice, species=struct.species, coords=struct.frac_coords, coords_are_cartesian=False)
        structures.append( (r, newStr))
    return structures