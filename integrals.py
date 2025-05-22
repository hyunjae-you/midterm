from scipy.special import erf
from scipy.special import gamma, gammainc
from hartree_fock.basis import get_h_basis
import numpy as np


def boys(n, x):
    """Boys 함수 계산"""
    if x < 1e-8:
        return 1.0 - x/3 + x**2/10  # 작은 x에 대한 Taylor 근사
    else:
        a = n + 0.5
        return gamma(a) * gammainc(a, x) / (2 * x**a)

def overlap_1D(xa, alpha, xb, beta):
    """1차원 가우시안 중첩 적분"""
    gamma_val = alpha + beta
    return np.sqrt(np.pi / gamma_val) * np.exp(-alpha*beta/gamma_val * (xa - xb)**2)

def overlap(a, b, Ra, Rb):
    """3차원 중첩 적분 (축별 계산)"""
    Sx = overlap_1D(Ra[0], a, Rb[0], b)
    Sy = overlap_1D(Ra[1], a, Rb[1], b)
    Sz = overlap_1D(Ra[2], a, Rb[2], b)
    return Sx * Sy * Sz

def nuclear_attraction(a, b, Ra, Rb, Rc, Z):
    p = a + b  # MATLAB의 p = a*b → Python에서 p = a + b로 수정
    Rp = (a * np.array(Ra) + b * np.array(Rb)) / p
    Rpc = np.linalg.norm(Rp - Rc)
    t = p * Rpc**2  # Boys 함수 입력값 일치화

    # 정규화 상수 반영 (MATLAB의 g1.N, g2.N)
    Na = (2*a/np.pi)**(3/4) * (4*a)**0  # s-orbital (l=0)
    Nb = (2*b/np.pi)**(3/4) * (4*b)**0
    V = -Z * 2*np.pi / p * boys(0, t) * overlap(a, b, Ra, Rb) * Na * Nb
    return V

def electron_repulsion(a, b, c, d, Ra, Rb, Rc, Rd, ca, cb, cc, cd):
    """4-중심 전자 반발 적분 (s-orbital 전용)"""
    # Gaussian 곱 정리 적용
    p = a + b
    pp = c + d
    Rp = (a * np.array(Ra) + b * np.array(Rb)) / p
    Rq = (c * np.array(Rc) + d * np.array(Rd)) / pp
    alpha_new = (p * pp) / (p + pp)
    Rpq = np.linalg.norm(Rp - Rq)
    T = alpha_new * Rpq**2

    # Boys 함수 계산
    F0 = boys(0, T)

    # 정규화 상수 및 수축 계수 반영
    Na = (2*a/np.pi)**(3/4) * (4*a)**0  # l=0 (s-orbital)
    Nb = (2*b/np.pi)**(3/4) * (4*b)**0
    Nc = (2*c/np.pi)**(3/4) * (4*c)**0
    Nd = (2*d/np.pi)**(3/4) * (4*d)**0
    norm = Na * Nb * Nc * Nd
    coeff = ca * cb * cc * cd

    # 전자 반발 적분 계산
    prefactor = (2 * np.pi**2.5) / (p * pp * np.sqrt(p + pp))
    integral = prefactor * F0 * overlap(a, b, Ra, Rb) * overlap(c, d, Rc, Rd)
    integral *= norm * coeff

    return integral

def compute_electron_repulsion(molecule, basis):
    num_atoms = len(molecule.atoms)
    num_basis = 4 * num_atoms
    eri = np.zeros((num_basis, num_basis, num_basis, num_basis))
    basis_set = get_h_basis()
    exponents = basis_set['exponents']
    coeffs = basis_set['coefficients']
    
    for p in range(num_basis):
        atom_p = p // 4
        Ra = molecule.atoms[atom_p][1]
        for q in range(num_basis):
            atom_q = q // 4
            Rb = molecule.atoms[atom_q][1]
            for r in range(num_basis):
                atom_r = r // 4
                Rc = molecule.atoms[atom_r][1]
                for s in range(num_basis):
                    atom_s = s // 4
                    Rd = molecule.atoms[atom_s][1]
                    integral = 0.0
                    for a in range(4):
                        alpha = exponents[a]
                        ca = coeffs[a]
                        for b in range(4):
                            beta = exponents[b]
                            cb = coeffs[b]
                            for c in range(4):
                                gamma = exponents[c]
                                cc = coeffs[c]
                                for d in range(4):
                                    delta = exponents[d]
                                    cd = coeffs[d]
                                    # 누락된 계수(ca, cb, cc, cd)를 인자에 추가
                                    integral += ca * cb * cc * cd * electron_repulsion(
                                        alpha, beta, gamma, delta,
                                        Ra, Rb, Rc, Rd,
                                        ca, cb, cc, cd  # 추가된 계수 인자
                                    )
                    eri[p, q, r, s] = integral
    return eri
def kinetic(a, b, Ra, Rb):
    """
    가우시안 기저 함수 간 운동 에너지 적분 계산 (s-orbital 전용)
    a, b: Gaussian 지수
    Ra, Rb: 원자 중심 좌표 (튜플)
    """
    # 가우시안 곱 정리 적용
    p = a + b
    Rp = (a * np.array(Ra) + b * np.array(Rb)) / p
    R_AB = np.linalg.norm(np.array(Ra) - np.array(Rb))
    q = (a * b) / p
    
    # 중첩 적분 계산
    S = overlap(a, b, Ra, Rb)
    
    # 운동 에너지 항 계산
    term1 = 3 * q * S
    term2 = 2 * q**2 * R_AB**2 * S
    T = term1 - term2
    
    return T

def compute_overlap(molecule, basis):
    num_atoms = len(molecule.atoms)
    num_basis = 4 * num_atoms  # 원자당 4개 기저 함수
    S = np.zeros((num_basis, num_basis))
    basis_set = get_h_basis()
    exponents = basis_set['exponents']
    coeffs = basis_set['coefficients']
    
    for i in range(num_basis):
        atom_i = i // 4
        Ra = molecule.atoms[atom_i][1]
        for j in range(num_basis):
            atom_j = j // 4
            Rb = molecule.atoms[atom_j][1]
            S_ij = 0.0
            for a in range(4):
                alpha = exponents[a]
                ca = coeffs[a]
                for b in range(4):
                    beta = exponents[b]
                    cb = coeffs[b]
                    S_ij += ca * cb * overlap(alpha, beta, Ra, Rb)
            S[i, j] = S_ij
    return S

def compute_kinetic(molecule, basis):
    num_atoms = len(molecule.atoms)
    num_basis = 4 * num_atoms
    T = np.zeros((num_basis, num_basis))
    basis_set = get_h_basis()
    exponents = basis_set['exponents']
    coeffs = basis_set['coefficients']
    
    for i in range(num_basis):
        atom_i = i // 4
        Ra = molecule.atoms[atom_i][1]
        for j in range(num_basis):
            atom_j = j // 4
            Rb = molecule.atoms[atom_j][1]
            T_ij = 0.0
            for a in range(4):
                alpha = exponents[a]
                ca = coeffs[a]
                for b in range(4):
                    beta = exponents[b]
                    cb = coeffs[b]
                    T_ij += ca * cb * kinetic(alpha, beta, Ra, Rb)
            T[i, j] = T_ij
    return T


def compute_nuclear(molecule, basis):
    num_atoms = len(molecule.atoms)
    num_basis = 4 * num_atoms
    V = np.zeros((num_basis, num_basis))
    basis_set = get_h_basis()
    exponents = basis_set['exponents']
    coeffs = basis_set['coefficients']
    
    for i in range(num_basis):
        atom_i = i // 4
        Ra = molecule.atoms[atom_i][1]
        for j in range(num_basis):
            atom_j = j // 4
            Rb = molecule.atoms[atom_j][1]
            V_ij = 0.0
            for a in range(4):
                alpha = exponents[a]
                ca = coeffs[a]
                for b in range(4):
                    beta = exponents[b]
                    cb = coeffs[b]
                    for (Z, Rc) in molecule.atoms:
                        integral = nuclear_attraction(alpha, beta, Ra, Rb, Rc, Z)
                        V_ij += ca * cb * integral
            V[i, j] = V_ij
    return V
