import atoms_shane as atoms
import vasp_shane as vasp
import numpy as np
import sys
def main():
    fn=sys.argv[1]
    ats = vasp.poscar_read(fn)
    vector=np.array(sys.argv[2].split(',')).astype(np.float)
    ats.translate_pos(vector)
    vasp.poscar_write(ats,"trans"+sys.argv[1])
main()
