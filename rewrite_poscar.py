import vasp_shane as vasp
import sys
fn = sys.argv[1]
g=vasp.poscar_read(fn)
g.return_ats_unit()
vasp.poscar_write(g,fn+".vasp")
