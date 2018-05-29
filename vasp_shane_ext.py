#this file includes useful routines based on the vasp_shane.py script. In this framework, atoms are saved as a dictionary for each element (as key) in the structure, and the fractional position of each atom as the value. Each should function by passing in an entire Atoms class, and return the Atoms class
import atoms_shane as atoms
import vasp_shane as vasp
import copy
def zero2one(Atoms): #simple script that returns any atoms outside of the unit cell to be between 0 and 1
    new_Atoms = copy.deepcopy(Atoms)
    new_Atoms.at = {x:[] for x in Atoms.at.keys()}
    for i in Atoms.at.keys():
        for j in Atoms.at[i]:
            new_Atoms.at[i].append([k % 1 for k in j])
    return new_Atoms
def remove_overlap(Atoms,tol): #remove first instance of overlapping atoms of the same species, within a tolerance specified in angstroms. Consider PBCs
