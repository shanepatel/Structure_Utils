import atoms_shane as atoms
import vasp_shane as vasp
import ascii_shane as ascii_s
import sys
f=vasp.poscar_read(sys.argv[1])
ascii_s.ascii_write(f,sys.argv[2])
