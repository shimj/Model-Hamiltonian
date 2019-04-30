import sympy as sp
from sympy.physics.quantum import TensorProduct
from itertools import combinations_with_replacement
from functools import reduce
import operator

pauli_matrices = [sp.eye(2), sp.Matrix([[0,1],[1,0]]), 
    sp.Matrix([[0,-sp.nsimplify("I")],[sp.nsimplify("I"),0]]), sp.Matrix([[1,0],[0,-1]])]

gamma_matrices = [TensorProduct(i,j) for i in pauli_matrices for j in pauli_matrices]
gamma_6m = [gamma_matrices[0], gamma_matrices[5], gamma_matrices[9], gamma_matrices[13], gamma_matrices[2], gamma_matrices[3]]

pg_matrices = [TensorProduct(i,j) for i in pauli_matrices for j in gamma_matrices]

def expr_num(n):
    assert n>=0
    return int((n+1)*n/2+n+1) # C_(n+1)^2 + C_(n+1)^1

def get_expr(n, symbolList=sp.symbols("kx ky kz")):
    assert n>=0
    if n==0: return [sp.nsimplify(1)]
    temp = combinations_with_replacement(symbolList, n)
    return [reduce(operator.mul, item) for item in temp]

def hermitian_inner_product(m1,m2):
    return (m1.H*m2).trace().simplify()

def get_hermitian_base(dim):
    """ get complete Hermitian basis with a certain dimension

    return (list of Matrix): a list of Hermitian basis

    dim (int): dimension of the Hermitian matrix
    """
    if dim==2:
        return pauli_matrices
    elif dim==4:
        return gamma_matrices
    elif dim==8:
        return pg_matrices
    assert dim > 0
    base = [sp.diag(*[0]*i,1,*[0]*(dim-i-1)) for i in range(dim)]
    for row in range(dim):
        for col in range(row+1, dim):
            basis = [[0 for i in range(dim)] for j in range(dim)]
            basis[row][col] = 1
            basis = sp.Matrix(basis)
            base.append(basis+basis.H)

            basis = [[0 for i in range(dim)] for j in range(dim)]
            basis[row][col] = sp.nsimplify("I")
            basis = sp.Matrix(basis)
            base.append(basis+basis.H)
    return base


def nullspace_inner(matrix, free_var, reduced, pivots, cols):
    """ function used by nullspace()"""
    vec = [0]*cols
    vec[free_var] = 1
    for piv_row, piv_col in enumerate(pivots):
        for pos in pivots[piv_row+1:] + (free_var,):
            vec[piv_col] -= reduced[piv_row, pos]
    return matrix._new(cols, 1, vec)

def mp_nullspace_inner(a):
    return nullspace_inner(*a)

def debug_print(*a, do_print=True):
    if do_print:
        print(*a)

def nullspace(matrix, simplify=sp.simplify, pool="", debug=False):
    """ solve b in Ab=0 (modified from the sympy nullspace)

    return (list of Matrix): solution base, element is column matrix

    matrix (Matrix): A in Ab=0
    simplify: passed to Matrix.rref() method
    pool (pool): if pool object passed in, then use it
    debug (boolean): if True, print the progress infomation
    """
    if pool and sp.__version__[0] == "0":
        print('Warning: have not yet added the ability of parallel calculation for nullspace() to the current version of sympy')
        pool = ''
    if not pool:
        debug_print("solving nullspace", do_print=debug)
        result = matrix.nullspace(simplify=sp.simplify)
    else:
        debug_print("solving rref", do_print=debug)
        reduced, pivots = matrix.rref(simplify=simplify)
        cols = matrix.cols
        free_vars = [i for i in range(matrix.cols) if i not in pivots]
        debug_print("solving nullspace", do_print=debug)
        result = pool.map(mp_nullspace_inner, [[matrix,
            free_var, reduced, pivots, cols] for free_var in free_vars])
        #result = [nullspace_inner(matrix, free_var, reduced, pivots, cols) for free_var in free_vars] 
    return result


def intersection(*bases, pool="", debug=False):
    """ get base of intersection space
    
    return (list of list of number): intersection base

    bases (list of list of list of number): several bases
    pool (pool): if pool object passed in, then use it
    debug (boolean): if True, print the progress infomation
    """
    if len(bases) == 0: return []
    if len(bases) == 1: return bases
    bases = sorted(bases, key=len)
    if len(bases) == 2:
        if len(bases[0]) == 0 or len(bases[1]) == 0: return []
        matrix_base0 = sp.Matrix(bases[0]).T
        matrix_base1 = sp.Matrix(bases[1]).T
        large_matrix = matrix_base0.row_join(-matrix_base1)
        intersection_base = [(matrix_base0*item[0:len(bases[0]),0]).T.tolist()[0] for item in nullspace(
            large_matrix, pool=pool)]
        #debug_print(len(intersection_base), len(intersection_base[0]), do_print=True)
        return intersection_base
    return intersection(intersection(*bases[0:2]), *bases[2:], pool=pool, debug=False)