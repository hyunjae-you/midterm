import numpy as np

class Molecule:
    def __init__(self):
        self.atoms = []  # (전하, [x, y, z]) 저장
        self.n_electrons = 0  # 총 전자 수 속성 추가
    
    def add_atom(self, charge, position):
        self.atoms.append( (charge, position) )
        self.n_electrons += charge  # 전하 합산으로 전자 수 계산
    
    def nuclear_repulsion(self):
        energy = 0.0
        for i in range(len(self.atoms)):
            Zi, Ri = self.atoms[i]
            for j in range(i+1, len(self.atoms)):
                Zj, Rj = self.atoms[j]
                R = np.linalg.norm(np.array(Ri) - np.array(Rj))
                energy += Zi * Zj / R
        return energy
