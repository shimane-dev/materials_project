#!/usr/bin/env python
from pymatgen.core import Structure
import numpy as np

FULL = "POSCAR_orth_2x2x1.vasp"
VAC  = "POSCAR_orth_2x2x1_LiVac1.vasp"

Sf = Structure.from_file(FULL)  # Li抜く前
Sv = Structure.from_file(VAC)   # Li抜いた後（= A）

# ---- 対策B：消えた Li（vacancy）を同定 ----
def frac_dist(p, q):
    d = np.array(p) - np.array(q)
    d -= np.round(d)  # 最近接像（分率）
    return np.linalg.norm(d)

li_f = [s.frac_coords for s in Sf if s.species_string == "Li"]
li_v = [s.frac_coords for s in Sv if s.species_string == "Li"]

missing = []
for p in li_f:
    if min(frac_dist(p, q) for q in li_v) > 5e-2:  # しきい値は0.05程度
        missing.append(p)

assert len(missing) == 1, "欠損Liが一意に特定できません（FULL/VACの対応を確認）"
x_vac, y_vac, z_vac = missing[0]

# ---- 対策A：vacancyに最も近い 3a(Me) を選ぶ（同じ層を優先）----
# Me候補を抽出（Li以外かつ O 以外 = 遷移金属）
me_ids = [i for i,s in enumerate(Sv) if s.species_string not in ("Li","O")]

# まず z の“層近さ”で絞り込む
def zdist(a, b):
    dz = (a - b) % 1.0
    dz = min(dz, 1.0 - dz)
    return dz
z_me = np.array([Sv[i].frac_coords[2] % 1.0 for i in me_ids])
z_v  = z_vac % 1.0
mask_same_layer = np.array([zdist(z, z_v) < 0.15 for z in z_me])  # 層幅は0.15くらい
cands = [me_ids[i] for i,m in enumerate(mask_same_layer) if m]
if not cands:
    # 同層で見つからなければ全候補から
    cands = me_ids

# 2D最近傍（a,b のみ）で最終決定
def frac_dist2d(p, q):
    d = np.array(p[:2]) - np.array(q[:2])
    d -= np.round(d)
    return np.linalg.norm(d)

target_me = min(cands, key=lambda i: frac_dist2d(Sv[i].frac_coords, [x_vac, y_vac, z_vac]))
assert frac_dist2d(Sv[target_me].frac_coords, [x_vac, y_vac, z_vac]) < 0.2, "対応3aが遠すぎます"

# ---- A/B を保存 ----
A = Sv.copy()  # A = Li欠損だけ入れた状態
B = Sv.copy()
species = Sv[target_me].species
B.remove_sites([target_me])
B.append(species, [x_vac, y_vac, z_vac], coords_are_cartesian=False)

A.to(fmt="poscar", filename="POSCAR_A")
B.to(fmt="poscar", filename="POSCAR_B")
print("Done: POSCAR_A / POSCAR_B")
