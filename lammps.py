import atoms_shane as atoms
import pickle
import os
dir_path = os.path.dirname(os.path.realpath(__file__))
def lammps_pos_write(atoms,filename,chg_dict={},atom_style='charges'):
    masses = pickle.load(open(dir_path+"/at_mass.p",'rb'))
    sat = sorted(list(atoms.at.keys()))
    nat = 0
    for i in sat:
        nat = nat + len(atoms.at[i])
    f = open(filename, "w")
    f.write("LAMMPS Description \n\n")
    f.write(str(nat) + " atoms\n\n")
    f.write(str(len(sat)) + " atom types\n\n")
    f.write("0  %10.8E  xlo xhi\n" %(atoms.cellvec[0][0]))
    f.write("0  %10.8E  ylo yhi\n" %(atoms.cellvec[1][1]))
    f.write("0  %10.8E  zlo zhi\n\n" %(atoms.cellvec[2][2]))
    f.write("Masses\n\n")
    for i in xrange(len(sat)):
        f.write(str(i+1)+" "+str(masses[sat[i]])+"\n")
    f.write("\nAtoms\n")
    ct = 1
    for i in xrange(len(sat)):
        for j in atoms.at[sat[i]]:
            if atom_style=='charges':
                f.write(str(ct)+"     "+str(i+1) + "       "+str(chg_dict[sat[i]]) + "   %9.5E   %9.5E   %9.5E\n" %(j[0],j[1],j[2]))
#            else:
#
            ct += 1
    f.close()
def lammps_pos_write_nonortho(atoms,filename,chg_dict={},atom_style='metal'):
    masses = pickle.load(open(dir_path+"/at_mass.p",'rb'))
    if atoms.coordinates=='Direct':
        atoms.swap_coordinate_system()
    sat = sorted(list(atoms.at.keys()))
    nat = 0
    for i in sat:
        nat = nat + len(atoms.at[i])
    f = open(filename, "w")
    xy = (atoms.cellvec[1][0]**2 + atoms.cellvec[0][1]**2)**0.5
    xz = (atoms.cellvec[2][0]**2 + atoms.cellvec[0][2]**2)**0.5
    yz = (atoms.cellvec[1][2]**2 + atoms.cellvec[2][1]**2)**0.5
    f.write("LAMMPS Description \n\n")
    f.write(str(nat) + " atoms\n\n")
    f.write(str(len(sat)) + " atom types\n\n")
    f.write("0  %10.8E  xlo xhi\n" %(atoms.cellvec[0][0]))
    f.write("0  %10.8E  ylo yhi\n" %(atoms.cellvec[1][1]))
    f.write("0  %10.8E  zlo zhi\n\n" %(atoms.cellvec[2][2]))
    f.write("%10.8E  %10.8E  %10.8E xy xz yz \n" %(atoms.cellvec[1][0]) %(atoms.cellvec[2][0]) %(atoms.cellvec[1][2]))
    f.write("Masses\n\n")
    for i in xrange(len(sat)):
        f.write(str(i+1)+" "+str(masses[sat[i]])+"\n")
    f.write("\nAtoms\n")
    ct = 1
    for i in xrange(len(sat)):
        for j in atoms.at[sat[i]]:
            if atom_style=='charges':
                f.write(str(ct)+"     "+str(i+1) + "       "+str(chg_dict[sat[i]]) + "   %9.5E   %9.5E   %9.5E\n" %(j[0],j[1],j[2]))
            elif atom_style=='metal':
                f.write(str(ct) + " " + str(ct) + " " + str(i+1) + " " + "0.0" + "   %9.5E   %9.5E   %9.5E\n" %(j[0],j[1],j[2]))
#
            ct += 1
    f.close()
def lammps_tschopp_read(filename):
    ats=atoms.Atoms()
    ats.coordinates="Cartesian"
    nat = 0
    rat = []
    sat = []
    Start = False
    Box = False
    Masses = False
    massvec = []
    specvec = []
    rat = []
    for line in f.readlines():
        if line.split() == []:
            continue
        if "atom types" in line:
            Start = True
            continue
        if Start == True and Box == False:
            if "xlo" in line:
                ats.cellvec[0][0] = float(line.split()[1])
                continue
            if "ylo" in line:
                ats.cellvec[1][1] = 2 * float(line.split()[1])
                hold = float(line.split()[1])
                continue
            if "zlo" in line:
                ats.cellvec[2][2] = float(line.split()[1])
                Box = True
                continue
        if Box == True: 
            if "Atoms" in line:
                Atoms = True
                continue
        if Masses == True:
            specvec.append(line.split()[1])
            rat.append([float(line.split()[-3]),float(line.split()[-2]),float(line.split()[-1])])
    sat = [masses[int(float(i))] for i in massvec]
    ats.at = {spec:[] for spec in sat}
    #    if line =="LAMMPS Description" or 
    for i in xrange(len(specvec)):
        ats.at[sat[int(specvec[i])-1]].append(rat[i])
    return ats

def lammps_pos_read(filename):
    f = open(filename,"r")
    masses = pickle.load(open(dir_path+"/mass_at.p",'rb'))
    ats=atoms.Atoms()
    ats.coordinates="Cartesian"
    nat = 0
    rat = []
    sat = []
    Start = False
    Box = False
    Masses = False
    massvec = []
    specvec = []
    rat = []
    for line in f.readlines():
        if line.split() == []:
            continue
        if "atom types" in line:
            Start = True
            continue
        if Start == True and Box == False:
            if "xlo" in line:
                ats.cellvec[0][0] = float(line.split()[1])
                continue
            if "ylo" in line:
                ats.cellvec[1][1] = float(line.split()[1])
                continue
            if "zlo" in line:
                ats.cellvec[2][2] = float(line.split()[1])
                Box = True
                continue
        if Box == True and Masses == False: 
            if "Masses" in line:
                continue
            if "Atoms" in line:
                Masses = True
                continue
            massvec.append(line.split()[1])
        if Masses == True:
            specvec.append(line.split()[1])
            rat.append([float(line.split()[-3]),float(line.split()[-2]),float(line.split()[-1])])
    sat = [masses[int(float(i))] for i in massvec]
    ats.at = {spec:[] for spec in sat}
    #    if line =="LAMMPS Description" or 
    for i in xrange(len(specvec)):
        ats.at[sat[int(specvec[i])-1]].append(rat[i])
    return ats
