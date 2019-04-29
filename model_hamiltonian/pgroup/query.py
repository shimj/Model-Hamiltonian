import sympy as sp
import copy
from itertools import combinations
from ._tools import accumulate, directSum

def chars_to_reps(Info, chars_including_E):
    if type(chars_including_E[0]) != type([0]):
        curr_dim = chars_including_E[0]
        curr_chars = chars_including_E[1:]
        for repsN in ["reps1", "reps2"]:
            jsonReps = Info[repsN]
            completed = False
            for repsi in jsonReps:
                if type(repsi["chars"][0]) == type([0]):
                    continue
                if repsi["d"]!=str(curr_dim):
                    continue
                ## simplify() can transform a string or a number into a sympy expression type.
                absDiffs = [abs(sp.simplify(ii).evalf()-sp.simplify(jj).evalf()) 
                    for ii, jj in zip(repsi["chars"], curr_chars)]
                if [ii<1e-8 for ii in absDiffs] == [True for ii in absDiffs]:
                    completed = True
                    if repsi["vecs_rep"]:
                        return copy.deepcopy(repsi["vecs_rep"])
                    else:
                        raise Exception("Data Missing")
        raise Exception("Characters Matching Failed")
    else:
        curr_dim = str(sum(sp.Matrix(chars_including_E[0])))
        curr_chars = list(map(list, zip(*chars_including_E[1:])))
        curr_chars_expand = [char for sublist in curr_chars for char in sublist]
        for repsN in ["reps1", "reps2"]:
            jsonReps = Info[repsN]
            completed = False
            for repsi in jsonReps:
                if type(repsi["chars"][0]) != type([0]):
                    continue
                if repsi["d"]!=str(curr_dim):
                    continue
                charsExpand = [char for sublist in repsi["chars"] for char in sublist]
                ## simplify() can transform a string or a number into a sympy expression type.
                absDiffs = [abs(sp.simplify(ii).evalf()-sp.simplify(jj).evalf()) 
                    for ii, jj in zip(charsExpand, curr_chars_expand)]
                if [ii<1e-8 for ii in absDiffs] == [True for ii in absDiffs]:
                    completed = True
                    if repsi["vecs_rep"]:
                        return copy.deepcopy(repsi["vecs_rep"])
                    else:
                        raise Exception("Data Missing")
        raise Exception("Characters Matching Failed")

def find_all_by_orbs_or_chars(Info, orbs=[], arrangement=[], chars_including_E=[]):
    """ find the items (one item contains not only a list of matrices) matching
    a single set of characters or orbitals or the arrangement, and then return
    all matched items. the data structures are a little different from that 
    of the original ones.
    In most cases, one should pass through either only characters or both 
    the orbitals and the arrangement. Note that the latter can only match one item.

    return (list of dict): all the matched item
    
    Info (dict): information about a certain group extracted from the json
    orbs (list of string): a list of atomic orbitals forming a irreducible representation
    arrangement (list of string): the characters of generators (excluding E) formed by
        the pure arrangement of atoms (i.e. ignoring the shape of the atomic orbitals)
    chars_including_E (list of string): the characters of generators and E, and the first
        one is the character of E i.e. the dimension of the irreducible representation.
    """
    if chars_including_E:
        vecs_reps1 = chars_to_reps(Info, chars_including_E)
    else:
        vecs_reps1 = []
        for repsN in ["reps1", "reps2"]:
            for reps in Info[repsN]:
                for item in reps["vecs_rep"]:
                    item = copy.deepcopy(item)
                    item["chars"] = reps["chars"]
                    vecs_reps1.append(item)
    if arrangement:
        arrangement = [str(item) for item in arrangement]
        vecs_reps2 = []
        for item in vecs_reps1:
            if sorted(item["arrangement"]) == sorted(arrangement):
                vecs_reps2.append(item)
    else:
        vecs_reps2 = vecs_reps1
    if orbs:
        vecs_reps3 = []
        for item in vecs_reps2:
            if sorted(item["vecs"]) == sorted(orbs):
                vecs_reps3.append(item)
    else:
        vecs_reps3 = vecs_reps2
    return vecs_reps3

def get_the_first_irred_rep(Info, orbs=[], arrangement=[], chars_including_E=[]):
    """ extract the irreducible representation (including the matrix of TR if exists) 
    from the first item among the data matching a single set of characters or 
    orbitals or the arrangement.
    In most cases, one should pass through either only characters or both 
    the orbitals and the arrangement. Note that the latter can only match one item.
    
    return (list of matrices): the first matched representation
    
    Info (dict): information about a certain group extracted from the json
    orbs (list of string): a list of atomic orbitals forming a irreducible representation
    arrangement (list of string): the characters of generators (excluding E) formed by
        the pure arrangement of atoms (i.e. ignoring the shape of the atomic orbitals)
    chars_including_E (list of string): the characters of generators and E, and the first
        one is the character of E i.e. the dimension of the irreducible representation.
    """
    data = find_all_by_orbs_or_chars(Info, orbs=orbs, 
        arrangement=arrangement, chars_including_E=chars_including_E)
    if len(data) == 0:
        return []
    else:
        if "repTR" in data[0]:
            return [*data[0]["rep"], data[0]["repTR"]]
        else:
            return data[0]["rep"]


def multi_orbs_to_rep(Info, multi_orbs_with_arr):
    multi_band_reps = []
    for one_band in multi_orbs_with_arr:
        multi_band_reps.append([sp.Matrix(item) for item 
            in get_the_first_irred_rep(Info, orbs=one_band[0], arrangement=one_band[1])])
    final_rep = [accumulate(directSum, item) for item in list(map(list, zip(*multi_band_reps)))]
    return final_rep

def multi_chars_to_rep(Info, chars_for_bands):
    multi_band_reps = []
    for chars_including_E in chars_for_bands:
        multi_band_reps.append([sp.Matrix(item) for item 
            in get_the_first_irred_rep(Info, chars_including_E=chars_including_E)])
    final_rep = [accumulate(directSum, item) for item in list(map(list, zip(*multi_band_reps)))]
    return final_rep

def match_info(Info, multi_orbs_with_arr=[], chars_for_bands=[], similar_trans=False):
    """ extract the irreducible representation (including the matrix of TR if exists) 
    from the first item among the data matching a single set of characters or 
    orbitals or the arrangement, and then pack them up with other important information.
    In most cases, one should pass through either only characters or both 
    the orbitals and the arrangement. Note that the latter can only match one item.
    
    return (list of dict): each dict describe a symmetry operation
        dict keys: "op3x3", "rep", "is_au", "name" (op3x3 is "gamma_nu" in the note, and it should be noted
            that the value of time reversal operator is a minus identity rather than the identity.)
    
    Info (dict): information about a certain group extracted from the json
    multi_orbs_with_arr (a list of list): each element describes one band (one irreducible representation)
        and is a two-element list: [orbs, arrangement] where orbs and arrangement have the same structure as
        defined in get_the_first_irred_rep()
        [0][0] (list of string): a list of atomic orbitals forming a irreducible representation
        [0][1] (list of string): the characters of generators (excluding E) formed by
            the pure arrangement of atoms (i.e. ignoring the shape of the atomic orbitals)
    chars_for_bands (list of string): each element describe one band and has the same structure as 
         the variable "chars_including_E" defined in get_the_first_irred_rep()
         [0] (list of string): the characters of generators and E, and the first
        one is the character of E i.e. the dimension of the irreducible representation.
    """
    if multi_orbs_with_arr:
        rep = multi_orbs_to_rep(Info, multi_orbs_with_arr=multi_orbs_with_arr)
    elif chars_for_bands:
        rep = multi_chars_to_rep(Info, chars_for_bands=chars_for_bands)
    else:
        raise Exception("Require at least one of multi_orbs_with_arr and chars_for_bands")
    if not rep: raise Exception("Matching failed")

    TR = True if "repTR" in str(Info) else False

    if similar_trans:
        if similar_trans.shape != rep[0].shape: raise Exception(
            "Wrong dimension of the similar-transform matrix")
        rep = [*[similar_trans.inv()*item*similar_trans for item in rep[0:-1]],
            similar_trans.inv()*rep[-1]*similar_trans.conjugate()] if TR else [
            similar_trans.inv()*item*similar_trans for item in rep]

    opers_name = [*Info["genes"], "TR"] if TR else Info["genes"]
    opers_3x3 = [*[sp.Matrix(item) for item in Info["3x3"]],
        -sp.eye(3)] if TR else [sp.Matrix(item) for item in Info["3x3"]]
    opers_is_au = [*[False]*len(Info["genes"]), True] if TR else [False]*len(Info["genes"])
    operations = []
    for item1, item2, item3, item4 in zip(opers_3x3, rep, opers_is_au, opers_name):
        operations.append({"op3x3":item1,"rep":item2,"is_au":item3,"name":item4})
    return operations



# other functionalities
def chars_to_orbs(Info, chars_including_E, arrangement=[]):
    data = find_all_by_orbs_or_chars(Info, chars_including_E=chars_including_E)
    result = []
    for item in data:
        if arrangement and sorted([str(item) for item in arrangement]) != sorted(item["arrangement"]):
            continue
        result.append([item["vecs"], item["arrangement"]])
    return result

def orbs_to_chars(Info, orbs, arrangement):
    data = find_all_by_orbs_or_chars(Info, orbs=orbs, arrangement=arrangement)
    result = []
    for item in data:
        if sorted([str(item) for item in arrangement]) == sorted(
            item["arrangement"]) and sorted(orbs) == sorted(item["vecs"]):
            return item["chars"]
    return []