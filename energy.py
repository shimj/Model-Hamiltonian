from model_hamiltonian._tool import get_hermitian_base
from model_hamiltonian.model_hamiltonian import energy
a = get_hermitian_base(4)
for n, item in enumerate(a):
    print(n, item)
b = [a[12], a[4], a[8]]

print(energy(b))