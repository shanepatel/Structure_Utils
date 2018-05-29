#!/usr/bin/env python
import sys
import atoms_shane
import copy
import vasp_shane
import ascii_shane

if len(sys.argv) < 3:
    print "usage: ascii2poscar.py input_filename output_filename"
    exit()
else:
    filename = sys.argv[1]
    writename = sys.argv[2]
atoms=ascii_shane.ascii_read(filename)
if atoms.coordinates=='Cartesian':
    atoms.swap_coordinate_system()
atoms.return_ats_unit()
vasp_shane.poscar_write(atoms,writename)
