#!/usr/bin/env python
from pymatgen.core import Structure, Site
import numpy as np

# 出発：Li層にLi欠損を1個入れた 2×2×1 構造（A）
A = Structure.from_file("POSCAR_orth_2x2x1_LiVac1.vasp")

# 1) Li vacancy の位置を推定
#   直前に Li を抜いたときの座標が分かるならそれを使うのが最良。
#   ここでは「Li層(z~0.5) で、近傍にLiが欠けている“穴”の幾何中心」を近似的にとる簡易法。
li_xyz = [site.frac_coords for site in A if site.species_string=="Li" and abs((site.frac_coords[2]%1)-0.5)<0.15]
# 面内グリッドで空きを推定するのは面倒なので、実用上は
#   → Li を抜いたときの (x_vac, y_vac, z_vac) をログ保存しておき、それを読み込むのが確実。
# ここでは例として x_vac, y_vac, z_vac を手で与える：
x_vac, y_vac, z_vac = 0.0, 0.0, 0.5   # ←実際の欠損座標に置き換えて！

# 2) vacancy と同じ in-plane の 3a(Me) を探す
def nearly(a,b,atol=1e-4):
    return np.allclose((np.array(a)-np.array(b))%1.0, [0,0], atol=atol) or \
           np.allclose((np.array(b)-np.array(a))%1.0, [0,0], atol=atol)

me_ids = [i for i,s in enumerate(A) if s.species_string!="Li" and abs(s.frac_coords[2]%1.0) < 0.08]  # 3a: z~0
target_me = None
for i in me_ids:
    xy = A[i].frac_coords[:2]
    if nearly(xy, [x_vac, y_vac], atol=2e-3):
        target_me = i
        break
assert target_me is not None, "対応する3a(Me)が見つかりません。x_vac,y_vacを確認してください。"

# 3) B 構造を作る：その Me を 3b(vac) へ“移動”（= 3a 側を空席化）
B = A.copy()
species = A[target_me].species
# 3a 側の Me を除去
B.remove_sites([target_me])
# vacancy 座標に Me を追加（結晶座標で追加）
B.append(species, [x_vac, y_vac, z_vac], coords_are_cartesian=False)

# 4) 保存
A.to(fmt="poscar", filename="POSCAR_A_3a_with_LiVac3b.vasp")       # 状態A
B.to(fmt="poscar", filename="POSCAR_B_Me_at_3b_vac_from_3a.vasp")  # 状態B
