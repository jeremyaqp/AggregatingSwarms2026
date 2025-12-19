"""
Converts .cfg files to .csv, while merging the (mux, muy) properties describing
the orientation vector into a single variable "theta".

"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys
from glob import glob


def cfg_to_csv(path_to_file, destination_folder):
    export_path = os.path.join(destination_folder, path_to_file.split("/")[-1][:-4] + ".csv")
    if os.path.exists(export_path): return "File exists : " + export_path
    with open(path_to_file, "r") as f:
        big_list = list()
        register = False
        part = 0
        columns = list()
        frame = 0
        for line in f.readlines():
            if "ATOMS id x y" in line:
                register = True
                columns = line.split()[2:]
                part = 0
                continue
            elif "ITEM" in line and register is True:
                register = False
                frame += 1
                continue
            if register:
                row = [float(p) for p in line.split()]
                if "mux" in columns:
                    theta = np.arctan2(row[columns.index("muy")], row[columns.index("mux")])
                    remove_elements = [row[columns.index("muy")], row[columns.index("mux")]]
                    row.remove(remove_elements[0])
                    row.remove(remove_elements[1])
                    row.append(theta)
                part = row[columns.index("id")]
                big_list.append(row + [frame, part])
    if "mux" in columns:
        columns.remove("mux")
        columns.remove("muy")
        df = pd.DataFrame(big_list, columns=columns + ["theta", "frame", "particle"])
    else:
        df = pd.DataFrame(big_list, columns=columns + ["frame", "particle"])
    df.to_csv(export_path, index=False)
    return export_path


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python cfg_to_csv.py <folder_path>")
        sys.exit(1)
    folder_dest = os.path.join(sys.argv[1], "CSV/")
    if not os.path.exists(folder_dest):
        os.makedirs(folder_dest)
    nb_done = len(os.listdir(folder_dest)) 
    nb_todo = len(glob(os.path.join(sys.argv[1], "*.cfg")))
    print(f"Converting {nb_todo} files, {nb_done} already done")
    if nb_done >= nb_todo:
        print(f"{sys.argv[1]} : All files already converted")
        sys.exit(0)
    for p in glob(os.path.join(sys.argv[1], "*.cfg")):
        p2 = cfg_to_csv(p, os.path.join(sys.argv[1], "CSV/"))
        print(p2)


