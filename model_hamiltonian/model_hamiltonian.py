import sympy as sp
from sympy.physics.quantum import TensorProduct
import numpy as np
from ._tool import hermitian_inner_product, get_hermitian_base, get_expr, expr_num, intersection, nullspace
import random
from functools import reduce
import operator

if sp.__version__[0] != "0":
    def linear_solve(A,b):
        '''solve linear equation Ax=b

        return (tuple of number): a tuple of column element of x
        A (list of list of number): A matrix
        b (list of number): b vector
        '''
        return list(sp.linsolve((sp.Matrix(A), sp.Matrix(b)), sp.symbols("x")))[0]
else:
    def linear_solve(A,b):
        '''solve linear equation Ax=b

        return (tuple of number): a tuple of column element of x
        A (list of list of number): A matrix
        b (list of number): b vector
        '''
        vnames = sp.symbols(["a"+str(i) for i in range(len(b))])
        result = sp.solve_linear_system(sp.Matrix(A).row_join(sp.Matrix(b)), *vnames)
        return tuple([result[item] for item in vnames])



def debug_print(*a, do_print=True):
    if do_print:
        print(*a)


#------------------------- matrix_M ----------------------------
def col_matrix_M(oper_rep, one_hermitian, hermitian_list, hermitian_list_norm2, is_au=False):
    '''calculate a certain column of matrix matrix_M

    return (list of number): a list of column element of matrix_M

    oper_op (Matrix): representation matrix of an operation (D(nv) or D(T))
    one_hermitian (Matrix): one of the Hermitian basis (B_i)
    hermitian_list (list of Matrix): a list of Hermitian basis
    hermitian_list_norm2 (list of number): a list of norm square of the Hermitian basis
    is_au (boolean): if oper_rep is antiunitary, you should pass true
    '''
    if is_au:
        new_hermitian = oper_rep * one_hermitian.conjugate() * oper_rep.H
    else:
        new_hermitian = oper_rep * one_hermitian * oper_rep.H
    return [hermitian_inner_product(new_hermitian, b) / norm2 for b, norm2 in zip(hermitian_list, hermitian_list_norm2)]
#print(col_matrix_M(sp.eye(2), pauli_matrices[0], pauli_matrices, [hermitian_inner_product(item, item) for item in pauli_matrices]))
def col_matrix_M_map(args):
    return col_matrix_M(*args)

def get_matrix_M(oper_rep, hermitian_list, is_au=False, pool=""):
    '''calculate the matrix matrix_M for a certain symmetry

    return (Matrix): the matrix matrix_M for the symmetry

    oper_op (Matrix): representation matrix of an operation (D(nv) or D(T))
    hermitian_list (list of Matrix): a list of Hermitian basis
    is_au (boolean): if oper_rep is antiunitary, you should pass true
    pool (pool): if pool object passed in, then use it
    '''
    hermitian_list_norm2 = [hermitian_inner_product(b, b) for b in hermitian_list]
    ## multiprocessing
    if pool:
        return sp.Matrix(pool.map(col_matrix_M_map, [[oper_rep, b, hermitian_list, hermitian_list_norm2, is_au] for b in hermitian_list])).transpose()
    else:
        return sp.Matrix([col_matrix_M(oper_rep, b, hermitian_list, hermitian_list_norm2, is_au) for b in hermitian_list]).transpose()


#--------------------------- matrix_F ----------------------------


def col_matrix_F(oper_3x3, oneExpr, expr_list, base_expr):
    '''calculate a certain column of matrix matrix_F

    return (list of number): a list of column element of matrix_F

    oper_3x3 (Matrix): 3x3 matrix of a point operation (tilde{g}_nu) or -diag(1,1,1) for TR
    oneExpr (Matrix): one of the expression basis (f_i(k))
    expr_list (list of Symbol): a list of expression basis (f_i(k))
    base_expr (list of Symbol): three symbols used in oneExpr and expr_list
    '''
    replace_list = list(zip(base_expr, sp.Matrix(oper_3x3).T * sp.Matrix(base_expr)))
    new_expr = oneExpr.subs(replace_list, simultaneous=True)
    A, b = [], []
    while sp.Matrix(A).rank() < len(expr_list):
        replace_list = [[item, random.randint(-100,100)] for item in base_expr]
        A.append([item.subs(replace_list) for item in expr_list])
        b.append(new_expr.subs(replace_list))
    return linear_solve(A, b)
    '''for Pythonista on iOS
    Ab = A
    for i in range(len(A)):
        Ab[i].append(b[i])
    vnames = sp.symbols(["a"+str(i) for i in range(len(A))])
    result = sp.solve_linear_system(sp.Matrix(Ab), *vnames)
    col = [result[item] for item in vnames]
    return col
    '''
    
'''test
klist = sp.symbols("kx ky kz")
kx, ky, kz = klist
print(col_matrix_F(-sp.eye(3), klist[1]**2*sp.sqrt(3)+klist[2]**2, [kx**2,ky**2,kz**2], klist))
'''
def col_matrix_F_map(args):
    return col_matrix_F(*args)

def get_matrix_F(oper_3x3, expr_list, base_expr, pool=""):
    '''calculate the matrix matrix_F for a certain symmetry

    return (Matrix): the matrix matrix_F

    oper_3x3 (Matrix): 3x3 matrix of a point operation (tilde{g}_nu) or -diag(1,1,1) for TR
    expr_list (list of Symbol): a list of expression basis (f_i(k))
    base_expr (list of Symbol): three symbols used in oneExpr and expr_list
    pool (pool): if pool object passed in, then use it
    '''
    if pool:
        return sp.Matrix(pool.map(col_matrix_F_map, [[oper_3x3, b, expr_list, base_expr] for b in expr_list])).transpose()
    else:
        return sp.Matrix([col_matrix_F(oper_3x3, b, expr_list, base_expr) for b in expr_list]).transpose()
'''test
klist = sp.symbols("kx ky kz")
print(get_matrix_F(-sp.eye(3), klist, klist))
'''


def solution_space_basis(oper, matrix_M, expr_list, base_expr, pool="", debug=False):
    '''solve (matrix_F otimes matrix_M)*alpha=alpha for a certain symmetry 

    ####return (list of Matrix): a list of vectors (single column Matrix object)
    return (list of list of number): a list of lists (each sublist represent a vector)

    oper (dict): a dictionary with four keys "op3x3", "rep", "is_au" and "name"
    matrix_M (Matrix): matrix_M matrix
    expr_list (list of Symbol): a list of expression basis (f_i(k))
    base_expr (list of Symbol): three symbols used in expr_list
    pool (pool): if pool object passed in, then use it
    debug (boolean): if True, print the progress infomation
    '''
    debug_print(oper["name"], "start", do_print=debug)
    matrix_F = get_matrix_F(oper["op3x3"], expr_list, base_expr, pool=pool)
    debug_print(oper["name"], "got matrix_F, calculating (FxM)*alpha=alpha", do_print=debug)
    matrix_FM = TensorProduct(matrix_F, matrix_M)
    temp = matrix_FM - sp.eye(matrix_FM.shape[0])
    # curr_basis = np.array(mat.nullspace(simplify=sp.nsimplify)).tolist()
    # solution_basis = [m.T.tolist()[0] for m in nullspace(mat, simplify=sp.nsimplify, pool=pool)]
    solution_basis = nullspace(temp, pool=pool)
    debug_print(oper["name"], "end", do_print=debug)
    return solution_basis

def narrow_hermitian_basis(PT_M_matrix, hermitian_list, pool=""):
    '''get a new set of Hermitian basis based on PT symmetry

    return (list of Matrix): a set of Hermitian basis

    one_or_two_operators (list of dict): a list of dict for PT operation, the dict have four keys "op3x3", "rep", "is_au" and "name"
    hermitian_list (list of Matrix): the original Hermitian basis
    pool (pool): if pool object passed in, then use it
    '''
    solution_basis = nullspace(PT_M_matrix-sp.eye(PT_M_matrix.shape[0]), simplify=sp.nsimplify, pool=pool)
    new_hermitian_list = [reduce(operator.add, [coef * hermitian for coef,
        hermitian in zip(solution, hermitian_list)]) for solution in solution_basis]
    return new_hermitian_list

def model_detail(operations, order_list, pool="", debug=False):
    '''get infomation of the result model Hamiltonian (core part)

    return (tuple):
        tuple[0] (list of Matrix): solution basis (A^{(l)})
        tuple[1] (list of Symbol): expression basis (f_i(k))
        tuple[2] (list of Matrix): the Hermitian basis (B_i)
        # tuple[3] (list of Symbol): a list of coefficients (c_l)

    operations (list of dict): a list of dict for all symmetries, the dict have four keys "op3x3", "rep", "is_au" and "name"
    order_list (list of int): a list of order of the polynomial of k
    pool (pool): if pool object passed in, then use it
    debug (boolean): if True, print the progress infomation
    '''
    assert order_list
    klist = sp.symbols("kx ky kz")
    repr_dim = operations[0]["rep"].shape[0]
    hermitian_list = get_hermitian_base(repr_dim)

    PT_contained = False
    PT_operation = [item for item in operations if item["is_au"]==True and item["op3x3"]==sp.eye(3)]
    if PT_operation:
        PT_contained = True
        PT_M_matrix = get_matrix_M(PT_operation[0]["rep"], hermitian_list, is_au=True, pool=pool)
    else:
        P_operation = [item for item in operations if item["is_au"]==False and item["op3x3"]==-sp.eye(3)]
        T_operation = [item for item in operations if item["is_au"]==True and item["op3x3"]==-sp.eye(3)]
        if P_operation and T_operation:
            PT_operation = [P_operation[0], T_operation[0]]
            PT_M_matrix = get_matrix_M(P_operation[0]["rep"]*T_operation[0]["rep"],
                hermitian_list, is_au=True, pool=pool)
        else:
            PT_M_matrix = ""
    if PT_M_matrix:
        debug_print("narrowing the number of Hermitian basis", do_print=debug)
        hermitian_list = narrow_hermitian_basis(PT_M_matrix, hermitian_list, pool=pool)
        debug_print("the number of the current Hermitian basis is", len(hermitian_list), do_print=debug)
        if PT_contained:
            ## delete PT
            operations = [item for item in operations if not (item["is_au"]==True and item["op3x3"]==sp.eye(3))]
        else:
            ## delete TR
            operations = [item for item in operations if not (item["is_au"]==True and item["op3x3"]==-sp.eye(3))]

    dim_matrix_M = len(hermitian_list)
    matrix_M_list = []
    debug_print("calculating matrix_M for all symmetries", do_print=debug)
    for oper in operations:
        matrix_M_list.append(get_matrix_M(oper["rep"], hermitian_list, is_au=oper["is_au"], pool=pool))
        debug_print("got matrix_M of ", oper["name"], do_print=debug)

    dim_matrix_F_list = [expr_num(order) for order in order_list]
    expr_list_all = []
    matrix_A_basis_all = []
    for now, order in enumerate(order_list):
        debug_print("\norder:", order, do_print=debug)
        expr_list = get_expr(order, klist)
        solution_bases = [solution_space_basis(oper, matrix_M_list[i], expr_list, klist,
            pool=pool, debug=debug) for i, oper in enumerate(operations)]
        debug_print("Got alpha bases for all symmetries", do_print=debug)
        ## merging list is much faster than row_join of matrices, so intersection() is designed to handle list of list
        solution_bases = [[m.T.tolist()[0] for m in base] for base in solution_bases]
        alpha_base = intersection(*solution_bases, pool=pool, debug=debug)
        debug_print("got intersection of solution spaces", do_print=debug)
        # alpha_base = intersection_basis(solution_bases)
        matrix_A_basis = [sp.Matrix(np.array(item).reshape(dim_matrix_F_list[now], dim_matrix_M)) for item in alpha_base]
        matrix_A_basis_temp = [sp.zeros(sum(dim_matrix_F_list[0:now]), dim_matrix_M).col_join(m) for m in matrix_A_basis]
        matrix_A_basis_expanded = [m.col_join(sp.zeros(sum(dim_matrix_F_list[now+1:]), dim_matrix_M)) for m in matrix_A_basis_temp]
        matrix_A_basis_all.extend(matrix_A_basis_expanded)
        expr_list_all.extend(expr_list)
    debug_print("\n", do_print=debug)

    return matrix_A_basis_all, expr_list_all, hermitian_list

def simple_calc(operations, order_list, pool="", debug=False):
    """ calculate the Hamiltonian
    
    return (list of tuple): [(h_1(k), B_1), (h_2(k), B_2), ...]

    order_list (list of int): a list of order of the polynomial of k
    pool (pool): if pool object passed in, then use it
    debug (boolean): if True, print the progress infomation
    """
    matrix_A_basis, expr_list, hermitian_list = model_detail(operations, order_list=order_list, pool=pool, debug=debug)
    coefficients = [sp.symbols("c"+str(i+1)) for i in range(len(matrix_A_basis))]

    if matrix_A_basis:
        ## sum on l
        matrix_A = reduce(operator.add, [c*A for c, A in zip(coefficients, matrix_A_basis)])
        ## sum on j
        hk = [reduce(operator.add, [coef*expr for expr, coef in zip(expr_list, col)]) for col in matrix_A.T.tolist()]
        ## return a list with index i (do not sum on i)
        result = []        
        for item1, item2 in zip(hk, hermitian_list):
            if item1 != 0:
                result.append((item1, item2))
    else: ## alpha (or A) has no solution
        result = []
    # result = [(num+1, item1, item2) for num, item1, item2 in zip(range(len(hk)), hk, hermitian_list) if item1 != 0]
    return result

def energy(hermitian_list):
    assert hermitian_list
    final_matrix = reduce(operator.add, [sp.symbols("h"+str(i+1))*m for i, m in enumerate(hermitian_list)])
    return final_matrix.eigenvals()

def eigenstates(hermitian_list):
    assert hermitian_list
    final_matrix = reduce(operator.add, [sp.symbols("h"+str(i+1))*m for i, m in enumerate(hermitian_list)])
    eigeninfo = final_matrix.eigenvects()
    # simplify the vectors
    eigeninfo = [(deg_item[0], deg_item[1], [sp.simplify(m) for m in deg_item[2]]) for deg_item in eigeninfo]
    return eigeninfo