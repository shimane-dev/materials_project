#!/usr/bin/env python

import os
from mp_api.client import MPRester
from pymatgen.core import Structure
# from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

API_KEY = 'fciIrSX7HleBYjvFNfVJlI8JimuTpJrd'
MP_ID = "mp-25373"
CIF_FILE = "LiMnO2_mp-25373.cif"

if os.path.exists(CIF_FILE):
    struct = Structure.from_file(CIF_FILE)
else:
    with MPRester(API_KEY) as mpr:
        struct = mpr.get_structure_by_material_id(MP_ID)
    struct.to(fmt="cif", filename=CIF_FILE)

# R-3m のままPOSCARで保存
struct.to(fmt="poscar", filename="POSCAR_R3m.vasp")
