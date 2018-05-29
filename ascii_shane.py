import sys
import atoms_shane as atoms
from poscar_aims2ascii_acf import latvec2dproj
from atoms_shane import *
import numpy as np
#*****************************************************************************************
def ascii_write(atoms,filename):
    if atoms.coordinates == "Direct":
        atoms.swap_coordinate_system()
    sat = []
    rat =[]
    vol = np.linalg.det(atoms.cellvec)
    for i in sorted(atoms.at.keys()):
        for j in atoms.at[i]:
            rat.append(j)
            sat.append(i)
    nat = len(rat)
    #dproj,cellvec,rxyz=latvec2dproj(atoms.cellvec,rat,len(rat))
    if filename=="screen":
        print "%d" % atoms.nat
        #print "%s" % atoms.boundcond
        print "%24.15E%24.15E%24.15E" % (atoms.cellvec[0][0],atoms.cellvec[1][0],atoms.cellvec[1][1])
        print "%24.15E%24.15E%24.15E" % (atoms.cellvec[2][0],atoms.cellvec[2][1],atoms.cellvec[2][2])
        for i in range(atoms.nat):
            x=atoms.rat[i][0]
            y=atoms.rat[i][1]
            z=atoms.rat[i][2]
            print "%24.15E%24.15E%24.15E%5s" % (x,y,z,atoms.sat[i])
    else:
        f= open(filename,"w")
        f.write(" %d" % nat + "     0.0000E+02   0.0000E+02          energy(eV) = 0.000E+03 enthalpy(ev)= 0.0000E+03 ucvol(ang3) = %24.15E  \n" %vol)
        f.write("%24.15E%24.15E%24.15E\n" % (atoms.cellvec[0][0],atoms.cellvec[1][0],atoms.cellvec[1][1]))
        f.write("%24.15E%24.15E%24.15E\n" % (atoms.cellvec[2][0],atoms.cellvec[2][1],atoms.cellvec[2][2]))
        for iat in xrange(len(rat)):
            x=float(rat[iat][0])
            y=float(rat[iat][1])
            z=float(rat[iat][2])
            s=sat[iat] 
            f.write(" %12.8f%12.8f%12.8f%5s\n" % (x,y,z,s))
        f.close()
#*****************************************************************************************
def ascii_read(filename):
    f=open(filename,"r")
    atoms=Atoms()
    atoms.boundcond="bulk"
    atoms.coordinates="Cartesian"
    iline=0
    iline_tot=0
    nconf=0
    nat = 0
    rat = []
    sat = []
    for line in f.readlines():
        iline_tot+=1
        tt=str(line).strip()
        if tt[0]=='#': continue
        #print tt
        iline+=1
        if iline==1:
            atoms.epot=float(line.split()[1])*27.211385
            #print atoms.epot
            pass
        elif iline==2:
            atoms.cellvec[0][0]=float(line.split()[0])
            atoms.cellvec[1][0]=float(line.split()[1])
            atoms.cellvec[1][1]=float(line.split()[2])
        elif iline==3:
            atoms.cellvec[2][0]=float(line.split()[0])
            atoms.cellvec[2][1]=float(line.split()[1])
            atoms.cellvec[2][2]=float(line.split()[2])
        else:
            nat+=1
            atoms.bemoved.append("TTT")
            rat.append([])
            icol=0
            #loop over the elemets, split by whitespace
            for i in line.split():
                #print i
                icol+=1
                if icol<4:
                    rat[-1].append(float(i))
                elif icol<5:
                    sat.append(i)
    species = list(set(sat))
    atoms.at = {spec:[] for spec in species}
    for i in xrange(len(sat)):
        atoms.at[sat[i]].append(rat[i])
    f.closed
    return atoms
#*****************************************************************************************
def ascii_write_fixat(atoms,filename,condition_vec):
    if atoms.coordinates == "Direct":
        atoms.swap_coordinate_system()
    sat = []
    rat =[]
    for i in sorted(atoms.at.keys()):
        for j in atoms.at[i]:
            rat.append(j)
            sat.append(i)
    nat = len(rat)
    #dproj,cellvec,rxyz=latvec2dproj(atoms.cellvec,rat,len(rat))
    if filename=="screen":
        print "%d" % atoms.nat
        #print "%s" % atoms.boundcond
        print "%24.15E%24.15E%24.15E" % (atoms.cellvec[0][0],atoms.cellvec[1][0],atoms.cellvec[1][1])
        print "%24.15E%24.15E%24.15E" % (atoms.cellvec[2][0],atoms.cellvec[2][1],atoms.cellvec[2][2])
        for i in range(atoms.nat):
            x=atoms.rat[i][0]
            y=atoms.rat[i][1]
            z=atoms.rat[i][2]
            print "%24.15E%24.15E%24.15E%5s" % (x,y,z,atoms.sat[i])
    else:
        f= open(filename,"w")
        f.write("%d\n" % nat)
        f.write("%24.15E%24.15E%24.15E\n" % (atoms.cellvec[0][0],atoms.cellvec[1][0],atoms.cellvec[1][1]))
        f.write("%24.15E%24.15E%24.15E\n" % (atoms.cellvec[2][0],atoms.cellvec[2][1],atoms.cellvec[2][2]))
        for iat in xrange(len(rat)):
            x=float(rat[iat][0])
            y=float(rat[iat][1])
            z=float(rat[iat][2])
            tf = str(condition_vec[iat])
            s=sat[iat] 
            f.write(" %12.8f%12.8f%12.8f%5s%5s\n" % (x,y,z,s,tf))
        f.close()
