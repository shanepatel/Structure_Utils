import numpy as np
import atoms_shane as atoms
import vasp_shane as vasp
def stack_struc(at_1,at_2):
    new_Atoms = atoms.Atoms()
    new_Atoms.coordinates = 'Cartesian'
    if at_1.coordinates == 'Direct':
        at_1.swap_coordinate_system()
    if at_2.coordinates == 'Direct':
        at_2.swap_coordinate_system()
    a = at_1.cellvec[0][0]
    b = at_1.cellvec[1][1]
    c = at_1.cellvec[2][2] +  at_2.cellvec[2][2]
    new_Atoms.cellvec=[[a,0,0],[0,b,0],[0,0,c]]
    new_Atoms.at = at_1.at
    for i in at_2.at.keys():
        #add_cords = []
        if i not in at_1.at.keys():
            new_Atoms.at[i] = []
        for j in at_2.at[i]:
            print j
            new_coord = np.array(j) + np.array(at_1.cellvec[2])
            new_Atoms.at[i].append(new_coord.tolist())
    return(new_Atoms)
