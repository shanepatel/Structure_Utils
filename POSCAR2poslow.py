import atoms_shane as vasp
import ascii_shane as asci
import vasp_shane as vasp
import glob
files = glob.glob("*.vasp")
n = len(files)
for i in xrange(1,n+1):
    string = 'poslow'+'{:05d}'.format(i)+'.ascii'
    print files[i-1],' , ' string
    f=vasp.poscar_read(files[i-1])
    asci.ascii_write(f,string)
