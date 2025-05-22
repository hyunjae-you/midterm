midterm project

-코드 중에서 ERI(전자 반발 적분)의 복잡성에서 4차원 배열 계산을 진행하기 때문에 계산과정이 오래걸리는 편입니다. 기다려주시면 감사하겠습니다...

-공지하였던 거리 범위 (0.1 ~ 2.0 angstrom)에서는 그래프가 확연히 보이지 않지만, 근거리에서 에너지를 과대해석하여 나타난 scaling의 문제임을 알 수 있었습니다. 1.0~3.0 angstrom으로 나타냈을 때 에너지에 대한 곡선이 reference의 유튜브 매트랩 강의에서와 유사한 형태를 가지고 있음을 확인하였습니다.

파일경로를 위와 같이 만들어주시고 py 파일들을 위와 같이 만들어주시면 됩니다.

\HF_project\

├─ setup.py

├─ hartree_fock\

│   ├─ __init__.py

│   ├─ basis.py

│   ├─ molecule.py

│   └─ scf.py

├─ examples\

│   └─ h2_energy.py      

└─ results\              ← 여기에 결과 자동 저장

-다음으로 "cd (HF_project의 파일경로)" 를 입력해주시고, "python .\examples\h2_energy.py" 를 해주시면 잠시후 계산이 시작됩니다.

- 그래프의 범위를 바꾸고 싶다면, h2_energy.py를 입력하신다음 19번째 줄에서 distances = np.linspace(0.1, 2.0, 20) 에서 수치를 변경시켜주면 됩니다. "linespace((그래프 시작지점),(그래프 끝지점), (간격))" 이 되겠습니다.

- DFT계산과 최대한 유사하게 구현하였기에 scf.py 에서 13번째 줄에 def scf(molecule, basis, max_iter=100, tol=1e-6): 을 보시면 'max_iter'은 최대 반복 회수, 'tol'은 수렴 허용 오차로 계산의 정확도와 관련된 수치도 구현했습니다. 이에 대해서도 변경가능하며 계산의 정확도 및 시간과 연관이 있습니다.


