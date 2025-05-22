import numpy as np
from scipy.linalg import eigh
from .integrals import compute_overlap, compute_kinetic, compute_nuclear, compute_electron_repulsion

def sort_eigs(eigvecs, eigvals):
    idx = np.argsort(eigvals)
    return eigvecs[:, idx], eigvals[idx]

def build_density(C, n_electrons):
    n_occ = (n_electrons + 1) // 2
    return 2 * C[:, :n_occ] @ C[:, :n_occ].T

def scf(molecule, basis, max_iter=100, tol=1e-6):
    S = compute_overlap(molecule, basis)
    
    # 중첩 행렬 양의 정부호 강제 (수치적 안정성 추가)
    min_eig = np.min(np.linalg.eigvalsh(S))
    if min_eig < 1e-8:
        S += np.eye(S.shape[0]) * (abs(min_eig) + 1e-8)
    
    T = compute_kinetic(molecule, basis)
    V = compute_nuclear(molecule, basis)
    eri = compute_electron_repulsion(molecule, basis)
    H_core = T + V

    # 초기 밀도 행렬
    eigvals, eigvecs = eigh(H_core, S)
    C, _ = sort_eigs(eigvecs, eigvals)
    D = build_density(C, molecule.n_electrons)

    # SCF 루프
    electronic_energy = []
    for _ in range(max_iter):
        J = np.einsum("pqrs,rs->pq", eri, D)
        K = np.einsum("prqs,rs->pq", eri, D)
        F = H_core + J - 0.5 * K

        eigvals, eigvecs = eigh(F, S)
        C, _ = sort_eigs(eigvecs, eigvals)
        D_new = build_density(C, molecule.n_electrons)

        energy = 0.5 * np.sum((H_core + F) * D)
        electronic_energy.append(energy)

        if np.linalg.norm(D_new - D) < tol:
            break
        D = D_new

    total_energy = electronic_energy[-1] + molecule.nuclear_repulsion()
    return total_energy
