#D3d+T: Bi2Se3 or MnBi2Te4_AFM
from model_hamiltonian.pgroup.get import get_data
from model_hamiltonian.pgroup import query
from model_hamiltonian.model_hamiltonian import simple_calc, energy
import multiprocessing
import sympy as sp
import time

try_parallel = True
gname = "D3d"
Info = get_data(TR=True)[gname]
print(Info["genes"])

bandchars = []
bandchars.append(["2", "1", "0", "2"]) # characters of identity and generators printed above for one band
bandchars.append(["2", "1", "0", "-2"])

## print possible orbitals
for num, item in enumerate(bandchars):
    print("band #"+str(num))
    print("\n".join([str(i) for i in query.chars_to_orbs(Info, chars_including_E=item)]))
    print("\n")

#choose the orbitals or use the first possible group of orbitals as default
'''
multi_orbs_with_arr = []
multi_orbs_with_arr.append([['p, 3/2, 1/2', 'p, 3/2, -1/2'], [1,-1,-1]])
multi_orbs_with_arr.append([['p, 3/2, 1/2', 'p, 3/2, -1/2'], [1,1,1]])
'''

order_list = [0, 1, 2]
similar_trans = sp.Matrix([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]]
)*sp.diag(1, sp.simplify("I"), 1, -sp.simplify("I"))
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
    print(item["name"], item["rep"])
print()

result = simple_calc(operations, order_list=order_list, pool=pool, debug=False)
if pool: pool.close()

print("result:\n","\n".join([str(item) for item in result]))
from functools import reduce
import operator
final_without_E = reduce(operator.add, [item[1]*item[2] for item in result if item[2] != sp.eye(item[2].shape[0])])
print("final without E:\n["+",\n".join([str(row) for row in final_without_E.tolist()])+"]")
print("(run time: "+str(round(time.time() - start_time, 1))+"s)")
try:
    print("\neigenvalues:", energy([item[2] for item in result]))
except:
    exit()
