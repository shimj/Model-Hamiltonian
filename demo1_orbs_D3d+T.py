#D3d+T: Bi2Se3 or MnBi2Te4_AFM
from model_hamiltonian.pgroup.get import get_data
from model_hamiltonian.pgroup import query
from model_hamiltonian.model_hamiltonian import simple_calc, show_result
import multiprocessing
import sympy as sp
import time

result_pattern = "i" #"i"
kz_zero = False ## True for xy surface
pretty_print = True

try_parallel = False
gname = "D3d"
Info = get_data(TR=True)[gname]
print(Info["genes"])

multi_orbs_with_arr = []
multi_orbs_with_arr.append([['p, 3/2, 1/2', 'p, 3/2, -1/2'], [1,-1,-1]])
multi_orbs_with_arr.append([['p, 3/2, 1/2', 'p, 3/2, -1/2'], [1,1,1]])

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

result = simple_calc(operations, order_list=order_list, output_index=result_pattern, pool=pool, debug=False)
if pool: pool.close()
print("(run time: "+str(round(time.time() - start_time, 1))+"s)")
show_result(result, result_pattern=result_pattern, kz_zero=kz_zero, pretty_print=pretty_print)