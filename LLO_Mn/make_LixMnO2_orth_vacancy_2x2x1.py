#!/usr/bin/env python
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.transformations.standard_transformations import SupercellTransformation
import numpy as np

##########
# 1) conventional cell（symprecは厳しめ）
s = Structure.from_file("LiMnO2_mp-25373.cif")
s_conv = SpacegroupAnalyzer(s, symprec=1e-5).get_conventional_standard_structure()

##########
# 2) 2×2×1
P = [[2,0,0],[0,2,0],[0,0,1]]
s_sup = SupercellTransformation(P).apply_transformation(s_conv)

##########
# 3) 指定層の Li(3b, z≈0.5) を1つ抜く
li_ids = [i for i,site in enumerate(s_sup) if site.species_string=="Li"]

def fold(z): return z % 1.0

# 3層のうち「中段」に最も近い層を選ぶ
layer_keys = {}
for i in li_ids:
    z = round(fold(s_sup[i].frac_coords[2]), 3)
    layer_keys.setdefault(z, []).append(i)

target_key = min(layer_keys, key=lambda z: abs(z-0.5))
s_sup.remove_sites([layer_keys[target_key][0]])

# 保存
s_sup.to(fmt="poscar", filename="POSCAR_orth_2x2x1_LiVac1.vasp")
