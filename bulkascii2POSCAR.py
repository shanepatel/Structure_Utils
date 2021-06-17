import atoms_shane as vasp
import ascii_shane as asci
import vasp_shane as vasp
import glob
import sys
if len(sys.argv) < 2:
    files = glob.glob("*poslow*.ascii")
    n = len(files)
    for i in xrange(1,n+1):
        string = 'poslow'+'{:05d}'.format(i)+'.vasp'
        print string, files[i-1]
        f=asci.ascii_read(files[i-1])
        vasp.poscar_write(f,string)
else:
    files = glob.glob(sys.argv[1])
    for i in files:
        print i
        string = i.replace(".ascii",".vasp")
        f=asci.ascii_read(i)
        vasp.poscar_write(f,string)
