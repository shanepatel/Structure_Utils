import numpy as np
import atoms_shane as atoms
import vasp_shane as vasp
import copy
def stack_struc(at_1_in,at_2_in,axis = 2):
    at_1=copy.deepcopy(at_1_in)
    at_2=copy.deepcopy(at_2_in)
    at_1.return_ats_unit()
    at_2.return_ats_unit()
    new_Atoms = atoms.Atoms()
    new_Atoms.coordinates = 'Cartesian'
    if at_1.coordinates == 'Direct':
        at_1.swap_coordinate_system()
    if at_2.coordinates == 'Direct':
        at_2.swap_coordinate_system()
    add_vec = [0,0,0]
    add_vec[axis] = 1
    a = np.array(at_1.cellvec[0]) + add_vec[0] * np.array(at_2.cellvec[0])
    b = np.array(at_1.cellvec[1]) + add_vec[1] * np.array(at_2.cellvec[1])
    c = np.array(at_1.cellvec[2]) + add_vec[2] * np.array(at_2.cellvec[2])
    new_Atoms.cellvec=[a,b,c]
    new_Atoms.at = at_1.at
    for i in at_2.at.keys():
        #add_cords = []
        if i not in at_1.at.keys():
            new_Atoms.at[i] = []
        for j in at_2.at[i]:
            new_coord = np.array(j) + np.array(at_1.cellvec[axis])
            new_Atoms.at[i].append(new_coord.tolist())
    #new_Atoms.swap_coordinate_system()
    return(new_Atoms)
    #vasp.poscar_write(new_Atoms,'compos_trash.vasp')
