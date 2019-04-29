def accumulate(function, list, *args_of_func):
    temp = list.pop(0)
    for item in list:
        temp = function(temp, item, *args_of_func)
    return temp

import sympy as sp
def directSum(m1,m2):
    return m1.row_join(sp.zeros(m1.shape[0],m2.shape[1])).col_join(sp.zeros(m2.shape[0],m1.shape[1]).row_join(m2))

def printMatrix(m):
    [print(item) for item in m.tolist()]


