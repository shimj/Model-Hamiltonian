#D3d+T: Bi2Se3 or MnBi2Te4_AFM
from model_hamiltonian.pgroup.get import get_data
from model_hamiltonian.pgroup import query
from model_hamiltonian.model_hamiltonian import simple_calc, energy
import multiprocessing
import sympy as sp
import time
from sympy import pprint

try_parallel = True
gname = "Oh"
Info = get_data(TR=True)[gname]
print(Info["genes"])

multi_orbs_with_arr = []
multi_orbs_with_arr.append([['s, 1/2', 's, -1/2'], [1,1,1]])
multi_orbs_with_arr.append([['p, 3/2, 3/2', 'p, 3/2, -3/2', 'p, 3/2, 1/2', 'p, 3/2, -1/2'], [1,1,1]])
multi_orbs_with_arr.append([['p, 1/2, 1/2', 'p, 1/2, -1/2'], [1,1,1]])

#multi_orbs_with_arr.append([['p, 3/2, 1/2', 'p, 3/2, -1/2'], [1,1,1]])

order_list = [1]
# similar_trans = sp.Matrix([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]]
#     )*sp.diag(1, sp.simplify("I"), 1, -sp.simplify("I"))
temp_m1 = sp.eye(8)


#---------******------cannot give corrent result but temp_m1[0,0] = -1 can.----------********------------
temp_m1[0,0] = 1
similar_trans = temp_m1
# temp_m2 = sp.eye(8)

# # temp_m1[4,4] = 0
# # temp_m1[5,5] = 0
# # temp_m1[4,5] = 1
# # temp_m1[5,4] = 1

# # temp_m2[3,3] = 0
# # temp_m2[5,5] = 0
# # temp_m2[3,5] = 1
# # temp_m2[5,3] = 1
# similar_trans = temp_m1*temp_m2
# similar_trans[0,0] = 1
# pprint(similar_trans)


#-------------------------------------------------------------------------------


if "bandchars" not in locals(): bandchars = []
if "multi_orbs_with_arr" not in locals(): multi_orbs_with_arr = []
if "similar_trans" not in locals(): similar_trans = False

start_time = time.time()
if try_parallel:
    cores = multiprocessing.cpu_count()
    try:
        pool = multiprocessing.Pool(processes=cores)
    except:
        pool = ""
else:
    pool = ""


operations = query.match_info(Info, multi_orbs_with_arr=multi_orbs_with_arr,
    chars_for_bands=bandchars, similar_trans=similar_trans)
for item in operations:
    print(item["name"],":")
    pprint(item["rep"])
print()

result = simple_calc(operations, order_list=order_list, pool=pool, debug=True)
if not result:
    print("Vanishing Halimtonian.")
    exit()
if pool: pool.close()

print("result:\n","\n".join([str(i+1)+" "+str(item) for i, item in enumerate(result)]))
from functools import reduce
import operator
final_without_E = reduce(operator.add,
    [item[0]*item[1] for item in result if item[1] != sp.eye(item[1].shape[0])],
    sp.zeros(result[0][1].shape[0]))
print("final without E:\n["+",\n".join([str(row) for row in final_without_E.tolist()])+"]")
print("(run time: "+str(round(time.time() - start_time, 1))+"s)")
try:
    print("\neigenvalues:", energy([item[1] for item in result]))
except:
    exit()
