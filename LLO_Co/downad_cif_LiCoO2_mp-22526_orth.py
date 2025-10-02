#!/usr/bin/env python

import os
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from mp_api.client import MPRester

API_KEY = "fciIrSX7HleBYjvFNfVJlI8JimuTpJrd"
MP_ID = "mp-22526"
CIF_FILE = "LiCoO2_mp-22526.cif"

if os.path.exists(CIF_FILE):
    struct = Structure.from_file(CIF_FILE)
else:
    with MPRester(API_KEY) as mpr:
        struct = mpr.get_structure_by_material_id(MP_ID)
    struct.to(fmt="cif", filename=CIF_FILE)

# conventional cell にする(直行系)
conv = SpacegroupAnalyzer(struct, symprec=1e-5).get_conventional_standard_structure()
conv.to(fmt="poscar", filename="POSCAR_orth.vasp")
