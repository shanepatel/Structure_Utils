import numpy as np
import atoms_shane as atoms
import vasp_shane as vasp
def swap_species_at_int(atoms, num,spec1,spec2,tol=0.125,ax=0):
    sat = sorted(atoms.at.keys())
    old_ats_1 = atoms.at[spec1]
    old_ats_2 = atoms.at[spec2]
    candidate_ats_1 = [i for i in old_ats_1 if (i[0] < tol) or (abs(1.0 - i[0]) < tol)]
    candidate_ats_2 = [i for i in old_ats_1 if (abs(0.5 - i[0]) < tol) or (abs(i[0] -0.5) < tol)] 
    base_ats_1 = [i for i in old_ats_1 if i not in candidate_ats_1 and i not in candidate_ats_2]
    candidate_ats_3 = [i for i in old_ats_2 if (i[0] < tol) or (abs(1.0 - i[0]) < tol)]
    candidate_ats_4 = [i for i in old_ats_2 if (abs(0.5 - i[0]) < tol) or (abs(i[0] -0.5) < tol)] 
    base_ats_2 = [i for i in old_ats_2 if i not in candidate_ats_3 and i not in candidate_ats_4]
    L1 = len(candidate_ats_1)
    L2 = len(candidate_ats_2)
    L3 = len(candidate_ats_3)
    L4 = len(candidate_ats_4)
    rands1 = list(reversed(sorted(np.random.choice(L1,num,replace=False))))
    rands2 = list(reversed(sorted(np.random.choice(L2,num,replace=False))))
    rands3 = list(reversed(sorted(np.random.choice(L3,num,replace=False))))
    rands4 = list(reversed(sorted(np.random.choice(L4,num,replace=False))))
    for i in xrange(num):
        at_1 = candidate_ats_1.pop(rands1[i])
        at_2 = candidate_ats_3.pop(rands3[i])
        candidate_ats_1.append(at_2)
        candidate_ats_3.append(at_1)
    for i in xrange(num):
        at_1 = candidate_ats_2.pop(rands2[i])
        at_2 = candidate_ats_4.pop(rands4[i])
        candidate_ats_2.append(at_2)
        candidate_ats_4.append(at_1)
    merged1 = base_ats_1 + candidate_ats_1 + candidate_ats_2
    merged2 = base_ats_2 + candidate_ats_3 + candidate_ats_4
    atoms.at[spec1] = merged1
    atoms.at[spec2] = merged2
    return atoms
def convert_species_at_int(atoms, num,spec1,spec2,tol=0.125,ax=0): #turn spec1 into spec2
    sat = sorted(atoms.at.keys())
    old_ats_1 = atoms.at[spec1]
    old_ats_2 = atoms.at[spec2]
    candidate_ats_1 = [i for i in old_ats_1 if (i[0] < tol) or (abs(1.0 - i[0]) < tol)]
    candidate_ats_2 = [i for i in old_ats_1 if (abs(0.5 - i[0]) < tol) or (abs(i[0] -0.5) < tol)] 
    base_ats_1 = [i for i in old_ats_1 if i not in candidate_ats_1 and i not in candidate_ats_2]
    candidate_ats_3 = [i for i in old_ats_2 if (i[0] < tol) or (abs(1.0 - i[0]) < tol)]
    candidate_ats_4 = [i for i in old_ats_2 if (abs(0.5 - i[0]) < tol) or (abs(i[0] -0.5) < tol)] 
    base_ats_2 = [i for i in old_ats_2 if i not in candidate_ats_3 and i not in candidate_ats_4]
    L1 = len(candidate_ats_1)
    L2 = len(candidate_ats_2)
    L3 = len(candidate_ats_3)
    L4 = len(candidate_ats_4)
    rands1 = list(reversed(sorted(np.random.choice(L1,num,replace=False))))
    rands2 = list(reversed(sorted(np.random.choice(L2,num,replace=False))))
    rands3 = list(reversed(sorted(np.random.choice(L3,num,replace=False))))
    rands4 = list(reversed(sorted(np.random.choice(L4,num,replace=False))))
    for i in xrange(num):
        at_1 = candidate_ats_1.pop(rands1[i])
        candidate_ats_3.append(at_1)
    for i in xrange(num):
        at_1 = candidate_ats_2.pop(rands2[i])
        candidate_ats_4.append(at_1)
    merged1 = base_ats_1 + candidate_ats_1 + candidate_ats_2
    merged2 = base_ats_2 + candidate_ats_3 + candidate_ats_4
    atoms.at[spec1] = merged1
    atoms.at[spec2] = merged2
    return atoms
def del_species_at_int(atoms, num,spec,tol=0.125,ax=0):
    sat = sorted(atoms.at.keys())
    old_ats = atoms.at[spec]
    candidate_ats_1 = [i for i in old_ats if (i[0] < tol) or (abs(1.0 - i[0]) < tol)]
    candidate_ats_2 = [i for i in old_ats if (abs(0.5 - i[0]) < tol) or (abs(i[0] -0.5) < tol)] 
    base_ats = [i for i in old_ats if i not in candidate_ats_1 and i not in candidate_ats_2]
    L1 = len(candidate_ats_1)
    L2 = len(candidate_ats_2)
    rands1 = reversed(sorted(np.random.choice(L1,num,replace=False)))
    rands2 = reversed(sorted(np.random.choice(L2,num,replace=False)))
    for i in rands1:
        del candidate_ats_1[i]
    for i in rands2:
        del candidate_ats_2[i]
    merged = base_ats + candidate_ats_1 + candidate_ats_2
    atoms.at[spec] = merged
    return atoms
