#!/usr/bin/env python
import os
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

# Materials Project API
from mp_api.client import MPRester

# ここに自分のAPIキーを設定（環境変数 MP_API_KEY でも可）
API_KEY = "YOUR_API_KEY"
MP_ID = "mp-25373"  # LiMnO2 (R-3m)

# ローカル保存ファイル名
LOCAL_CIF = "LiMnO2_mp-25373.cif"

def get_structure():
    """CIFファイルがあればそれを読む。なければMPから取得して保存。"""
    if os.path.exists(LOCAL_CIF):
        print(f"[INFO] {LOCAL_CIF} を読み込みます")
        struct = Structure.from_file(LOCAL_CIF)
    else:
        print(f"[INFO] {LOCAL_CIF} が無いのでAPIから取得します")
        with MPRester(API_KEY) as mpr:
            struct = mpr.get_structure_by_material_id(MP_ID)
        struct.to(fmt="cif", filename=LOCAL_CIF)
        print(f"[INFO] CIF を保存しました → {LOCAL_CIF}")
    return struct

def save_poscars(struct):
    """R-3mと直交セルのPOSCARを保存"""
    # そのまま（R-3m）
    struct.to(fmt="poscar", filename="POSCAR_R3m")
    print("[INFO] 保存しました → POSCAR_R3m")

    # 直交セル（標準直交セル）
    sga = SpacegroupAnalyzer(struct, symprec=1e-5)
    conv = sga.get_conventional_standard_structure()
    conv.to(fmt="poscar", filename="POSCAR_ortho")
    print("[INFO] 保存しました → POSCAR_ortho")

if __name__ == "__main__":
    structure = get_structure()
    save_poscars(structure)