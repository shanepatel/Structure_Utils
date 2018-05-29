import atoms_shane as atoms
import vasp_shane as vasp
f=vasp.poscar_read('STO_GB.vasp')
f.swap_axis(0,2)
vasp.poscar_write(f,'STO_GB_swap.vasp')
