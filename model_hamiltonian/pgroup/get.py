import json
import copy
import os.path
dir_path = os.path.abspath(os.path.dirname(__file__))

def get_data(TR):
    if TR:
        file_path = os.path.join(dir_path, "json/GroupInfo_with_TR_full.json")
    else:
        file_path = os.path.join(dir_path, "json/GroupInfo_no_TR_full.json")
    with open(file_path) as f:
        data = json.loads(f.read())
    data_copy = copy.deepcopy(data)
    for key, detail in data_copy.items():
        if detail["name2"]:
            data[detail["name2"]] = data[key]
    return data