import vasp
import numpy as np
import pickle as pkl
import glob
def main(filein,system='STO',MHM=True):
    pots = pkl.load(open('potdict.p', 'rb'))
    atoms = vasp.poscar_read(filein)
    fileout = filein.replace('.vasp','.in')
    cellvec = '  '+str(atoms.cellvec[0][0])+' '+str(atoms.cellvec[1][1])+' '+str(atoms.cellvec[2][2])+' 90 90 90\n'
    f = open(fileout,'w')
    if MHM == True and system =='STO':
        f.write('socket opti dist conp\n gtol 0.0005\n gmax opt 0.0000005\n maxcyc 10000000\nsocket host\nsto_XXX/\nsocket inet 0\n')
    else if MHM = False:
        f.write('single\n')
    if system == 'STO':
        f.write('cutp 20.0 polynomial 18.0\n  cell\n')
    f.write(cellvec)
    f.write('  cartesian\n')
    for i in xrange(len(atoms.rat)):
        linetext = str(atoms.sat[i]) +' '+str(atoms.rat[i][0])+' '+str(atoms.rat[i][1])+' '+str(atoms.rat[i][2])+'\n'
        f.write(linetext)
    f.write('\n'+pots['mccoy']+'\n')
    f.close()
poscars=glob.glob('*.vasp')
main('POSCAR.vasp')
