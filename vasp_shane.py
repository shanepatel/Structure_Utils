import atoms_shane as atoms
import os 
import copy
import subprocess
import re
#*****************************************************************************************
def poscar_write(atoms,filename):
    sat = sorted(atoms.at.keys())
    str_line = " ".join(sat)
    if filename=="screen":
        print str_line
        print "1.0"
        print "%24.15E%24.15E%24.15E" % (atoms.cellvec[0][0],atoms.cellvec[0][1],atoms.cellvec[0][2])
        print "%24.15E%24.15E%24.15E" % (atoms.cellvec[1][0],atoms.cellvec[1][1],atoms.cellvec[1][2])
        print "%24.15E%24.15E%24.15E" % (atoms.cellvec[2][0],atoms.cellvec[2][1],atoms.cellvec[2][2])
        print str_line
        str_line=""
        for i in sat:
            str_line+=" %d" % len(atoms.at[i])
        print str_line
        print atoms.coordinates
        for s in sat:
            for iat in atoms.at[s]:
                x=float(iat[0])
                y=float(iat[1])
                z=float(iat[2])
                print "%24.15E%24.15E%24.15E" % (x,y,z)
    else:
        f=open(filename,"w")
        f.write("%s\n" % str_line)
        f.write("1.0\n")
        f.write("{:.6f}  {:.6f}  {:.6f}".format(atoms.cellvec[0][0],atoms.cellvec[0][1],atoms.cellvec[0][2])+"\n")
        f.write("{:.6f}  {:.6f}  {:.6f}".format(atoms.cellvec[1][0],atoms.cellvec[1][1],atoms.cellvec[1][2])+"\n")
        f.write("{:.6f}  {:.6f}  {:.6f}".format(atoms.cellvec[2][0],atoms.cellvec[2][1],atoms.cellvec[2][2])+"\n")
        f.write("%s\n" % str_line)
        str_line=""
        for i in sat:
            str_line+=" %d" % len(atoms.at[i])
        f.write("%s\n" % str_line)
        f.write("%s\n" % atoms.coordinates)
        for s in sat:
            for iat in atoms.at[s]:
                x=float(iat[0])
                y=float(iat[1])
                z=float(iat[2])
                f.write("   {:.6f}  {:.6f}  {:.6f}".format(x,y,z)+"\n")
        f.close()
#*****************************************************************************************
def xdatcar_read(coordinates):
    f=open("XDATCAR","r")
    atoms_all=[]
    iline=0
    nconf=0
    for line in f.readlines():
        iline+=1
        if iline==1:
            atoms=atoms.Atoms()
            comment1=line
        elif iline==2:
            scale=float(line)
        elif iline==3:
            atoms.cellvec[0][0]=float(line.split()[0])
            atoms.cellvec[0][1]=float(line.split()[1])
            atoms.cellvec[0][2]=float(line.split()[2])
        elif iline==4:
            atoms.cellvec[1][0]=float(line.split()[0])
            atoms.cellvec[1][1]=float(line.split()[1])
            atoms.cellvec[1][2]=float(line.split()[2])
        elif iline==5:
            atoms.cellvec[2][0]=float(line.split()[0])
            atoms.cellvec[2][1]=float(line.split()[1])
            atoms.cellvec[2][2]=float(line.split()[2])
        elif iline==6:
            ntypes=0
            atom_types=[]
            for i in line.split():
                atom_types.append(i)
                ntypes+=1
        elif iline==7:
            ncol=0
            atoms.nat=0
            atom_nn=[]
            for i in line.split():
                atom_nn.append(int(i))
                atoms.nat+=int(i)
                ncol+=1
            if ncol!=ntypes:
                print "ERROR: in lines 6 and 7, number of columns differ: %3d%3d",ntypes,ncol
            for itype in range(ntypes):
                for iat in range(atom_nn[itype]):
                    atoms.sat.append(atom_types[itype])
        elif iline>7:
            if (iline-8)%(atoms.nat+1)>=1:
                atoms.rat.append([])
                xred=float(line.split()[0])
                yred=float(line.split()[1])
                zred=float(line.split()[2])
                if(coordinates=="Cartesian"):
                    atoms.rat[-1].append(xred*atoms.cellvec[0][0]+yred*atoms.cellvec[1][0]+zred*atoms.cellvec[2][0])
                    atoms.rat[-1].append(xred*atoms.cellvec[0][1]+yred*atoms.cellvec[1][1]+zred*atoms.cellvec[2][1])
                    atoms.rat[-1].append(xred*atoms.cellvec[0][2]+yred*atoms.cellvec[1][2]+zred*atoms.cellvec[2][2])
                elif(coordinates=="Reduced"):
                    atoms.rat[-1].append(xred)
                    atoms.rat[-1].append(yred)
                    atoms.rat[-1].append(zred)
                else:
                    print "ERROR: unknown coordinates"
                atoms.coordinates=coordinates
            if (iline-7)%(atoms.nat+1)==0:
                atoms_all.append(atoms.Atoms())
                atoms_all[-1]=copy.copy(atoms)
                nconf+=1
                atoms.rat=[]
    f.closed
    return atoms_all
#*****************************************************************************************
def poscar_read(filename):
    f=open(filename,"r")
    sat = []
    rat = []
    iline=0
    for line in f.readlines():
        iline+=1
        if iline==1:
            ats=atoms.Atoms()
            #nconf+=1
            comment1=line
        elif iline==2:
            comment2=line
        elif iline==3:
            ats.cellvec[0][0]=float(line.split()[0])
            ats.cellvec[0][1]=float(line.split()[1])
            ats.cellvec[0][2]=float(line.split()[2])
        elif iline==4:
            ats.cellvec[1][0]=float(line.split()[0])
            ats.cellvec[1][1]=float(line.split()[1])
            ats.cellvec[1][2]=float(line.split()[2])
        elif iline==5:
            ats.cellvec[2][0]=float(line.split()[0])
            ats.cellvec[2][1]=float(line.split()[1])
            ats.cellvec[2][2]=float(line.split()[2])
        elif iline==6:
            ntypes=0
            atom_types=[]
            for i in line.split():
                atom_types.append(i)
                ntypes+=1
        elif iline==7:
            ncol=0
            ats.nat=0
            atom_nn=[]
            for i in line.split():
                atom_nn.append(int(i))
                ats.nat+=int(i)
                ncol+=1
            if ncol!=ntypes:
                print "ERROR: in lines 6 and 7, number of columns differ: %3d%3d",ntypes,ncol
            for itype in range(ntypes):
                for iat in range(atom_nn[itype]):
                    sat.append(atom_types[itype])
        elif iline==8:
            if "direct" in line.strip().lower():
                ats.coordinates="Direct"
            else:
                ats.coordinates="Cartesian"
        elif iline>8:
            rat.append([])
            xred=float(line.split()[0])
            yred=float(line.split()[1])
            zred=float(line.split()[2])
            if(ats.coordinates=="Direct"):
                rat[-1].append(xred)
                rat[-1].append(yred)
                rat[-1].append(zred)
            elif(ats.coordinates=="Cartesian"):
                rat[-1].append(xred)
                rat[-1].append(yred)
                rat[-1].append(zred)
            else:
                print "ERROR: unknown coordinates"
            if not line.strip(): break
    species = list(set(sat))
    ats.at = {spec:[] for spec in species}
    for i in xrange(len(sat)):
        ats.at[sat[i]].append(rat[i])
    f.close()
    return ats
def read_outcar_energy(filename='OUTCAR',last=True):
    #command = "grep TOTEN " + filename + "| tail -1 | awk '{print $5}'"
    out= open(filename,"r")
    TOTENS=[]
    for line in out:
        if "TOTEN" in line:
            e = line.split(" ")[-2]
            TOTENS.append(e)
    #e = subprocess.Popen(command,shell=True)
    #e=os.system(command)
    if last==True:
        return(float(TOTENS[-1]))
#*****************************************************************************************
