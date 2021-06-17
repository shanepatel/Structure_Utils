import numpy as np
import atoms_shane as atoms
import vasp_shane as vasp
def add_vac(at,vac_thickness,axis = 2, opt = 2):#default is to add to c, pad symetrically. Vac thickness in angstroms
    if at.coordinates == 'Direct':
        at.swap_coordinate_system()
    if opt==2: #symmetric vacuum
    at.cellvec[2][2] = at.cellvec[2][2] + 2 * vac_thickness
    for i in at.at.keys():
        at.at[i] = at.at[i]+thickness
    return at
