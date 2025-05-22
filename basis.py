import numpy as np

def get_h_basis():
    """STO-3G 기저 세트 (수소 원자 전용)"""
    sto3g = {
        'exponents': [18.7311370, 2.8253937, 0.6401217, 0.1612778],
        'coefficients': [0.0334946, 0.2347269, 0.81375733, 1.0]
    }
    l = 0  # s-orbital
    norm = [(2 * ex / np.pi)**(3/4) for ex in sto3g['exponents']]
    sto3g['coefficients'] = [c * n for c, n in zip(sto3g['coefficients'], norm)]
    return sto3g
