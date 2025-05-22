import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import matplotlib.pyplot as plt
from hartree_fock.molecule import Molecule
from hartree_fock.scf import scf
from hartree_fock.basis import get_h_basis

basis = get_h_basis()  # 수소 기저 세트 불러오기

# 단일 수소 원자 에너지 계산
h_atom = Molecule()
h_atom.add_atom(1, (0,0,0))
h_energy = scf(h_atom, basis)

# H₂ 분자 결합 에너지 계산
distances = np.linspace(0.1, 2.0, 20)
binding_energies = []

print("\n===== H₂ 결합 에너지 계산 시작 =====")
for idx, r in enumerate(distances):
    h2 = Molecule()
    h2.add_atom(1, (0, 0, 0))
    h2.add_atom(1, (r, 0, 0))
    energy = scf(h2, basis)
    binding_energy = energy - 2 * h_energy
    binding_energies.append(binding_energy)
    print(f"진행률: {idx+1}/{len(distances)} | 거리: {r:.2f} Å | 결합 에너지: {binding_energy:.6f} Ha")

# h_energy 계산 코드 (examples/h2_energy.py)
h_atom = Molecule()
h_atom.add_atom(1, (0,0,0))
print(f"[검증] H 원자의 전자 수: {h_atom.n_electrons}")  # 1 출력 확인
h_energy = scf(h_atom, basis)  # -0.5 Hartree 정도의 값이 나와야 함
print(f"[결과] H 원자 에너지: {h_energy:.6f} Ha")  # ≈ -0.5 Ha 확인
print(f"H 원자 에너지: {h_energy}")  # 값 확인

# 그래프 그리기
plt.plot(distances, binding_energies)
plt.xlabel('Distance (Å)')
plt.ylabel('Binding Energy (Ha)')
plt.savefig('results/h2_plot.png')




