import atoms_shane as vasp
import ascii_shane as asci
import vasp_shane as vasp
import glob
files = glob.glob("*poslow*.ascii")
n = len(files)
for i in xrange(1,n+1):
    string = 'poslow'+'{:05d}'.format(i)+'.vasp'
    print string, files[i-1]
    f=asci.ascii_read(files[i-1])
    vasp.poscar_write(f,string)
