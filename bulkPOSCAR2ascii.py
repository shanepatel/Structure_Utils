import atoms_shane as vasp
import ascii_shane as asci
import vasp_shane as vasp
import glob
files = glob.glob("*.vasp")
n = len(files)
for i in files:
    string=i.replace('.vasp','.ascii')
    f=vasp.poscar_read(i)
    asci.ascii_write(f,string)
