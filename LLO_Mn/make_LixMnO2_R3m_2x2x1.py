#!/usr/bin/env python
from pymatgen.core import Structure
from pymatgen.transformations.standard_transformations import SupercellTransformation

s_conv = Structure.from_file("POSCAR_R3m.vasp")

P = [[2,0,0],[0,2,0],[0,0,1]]
#P = [[1,0,0],[0,2,0],[0,0,2]]
s_sup = SupercellTransformation(P).apply_transformation(s_conv)

s_sup.to(fmt="poscar", filename="POSCAR_R3m_2x2x1.vasp")
# s_sup.to(fmt="poscar", filename="POSCAR_R3m_1x2x2.vasp")
