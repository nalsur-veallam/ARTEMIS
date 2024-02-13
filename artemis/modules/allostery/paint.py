import sys
import json
import numpy as np
from pymol import cmd
from tqdm import tqdm
import pandas as pd

visual = "Y"
act = True
alls = True
noseq = 0
k = 0

def max_top(array, top):
    size = NResidues
    top10 = int(top/100*size)

    supp = np.zeros(size)
    for i in range(top10):
        supp[i] = 1

    d = {'data':array, 'index':np.arange(0,NResidues)}

    df = pd.DataFrame(data=d)
    df = df.sort_values(by=['data'], ascending=False)

    for i in range(NResidues):
        if np.array(df['data'])[i] < 0:
            supp[i] = 1

    df['data'] = df['data']*supp
    df = df.sort_values(by=['index'])
    return np.array(df['data'])

def draw_allostery(Allostery, out_path, top, noseq, str_path):

    print("\nSCRIPT FOR DISPLAYING ALLOSTERIC COMMUNICATION INTENSITY ON THE STRUCTURE IS LAUNCHED\n")

    if Allostery.allosteric_site is None:
        alls = False
    else:
        alls_arr = []
        group_alls = []
        group_alls_names = []

        for i in range(NResidues):
            if i + 1 in Allostery.allosteric_site:
                alls_arr.append(1)
            else:
                alls_arr.append(0)

    if top is None:
        intensity = []
        for i in range(Allostery.NResidues):
            if i + 1 in Allostery.active_site:
                intensity.append(-1)
            else:
                inten = 0
                for resid in Allostery.active_site:
                    if np.abs(resid - 1 - i) >= noseq:
                        inten += Allostery.map_[resid - 1][i]
                intensity.append(inten)
    else:
        intensity = np.zeros(Allostery.NResidues)

        for resid in Allostery.active_site:
            inten = []
            for i in range(Allostery.NResidues):
                if i + 1 in Allostery.active_site:
                    inten.append(-1)
                elif np.abs(resid - 1 - i) >= noseq:
                    inten.append(Allostery.map_[resid - 1][i])
                else:
                    inten.append(0)
            intensity += max_top(inten, top)


    intensity = np.array(intensity)
    intensity[intensity < 0] = -1

    resids = []
    group = []
    group_names = []

    def renumbering(resi, resn, b):
        resid = resi + resn
        if not(resid in resids):
            resids.append(resid)
            global k
            k += 1

            if intensity[k-1] < 0:
                group.append(resi)
                group_names.append(resn)

            if alls:
                if alls_arr[k-1] == 1:
                    group_alls.append(resi)
                    group_alls_names.append(resn)

            return intensity[k-1]
        return intensity[k-1]

    myspace = {'renumbering': renumbering}

    cmd.load(str_path)
    cmd.alter("all",'b = renumbering(resi, resn, b)', space=myspace)

    if visual=="Y":
        obj = cmd.get_object_list("all")[0]
        cmd.show_as("cartoon",obj)
        cmd.cartoon("putty", obj)
        cmd.set("cartoon_putty_scale_min", min(intensity),obj)
        cmd.set("cartoon_putty_scale_max", max(intensity),obj)
        cmd.set("cartoon_putty_transform", 0,obj)
        cmd.set("cartoon_putty_radius", 0.4,obj)
        cmd.spectrum("b","rainbow", obj)
        cmd.ramp_new("count", obj, [min(intensity), max(intensity)], "rainbow")
        cmd.recolor()

    for j in range(len(group)):
        text = "resi " + str(group[j]) + "& resn " + str(group_names[j])
        cmd.select("active_site", text, merge=1)

    if act:
        obj = "active_site"
        cmd.show_as("spheres",obj)
        cmd.color("white", obj)
        cmd.recolor()

    if alls:
        for j in range(len(group_alls)):
            text = "resi " + str(group_alls[j]) + "& resn " + str(group_alls_names[j])
            cmd.select("allosteric_site", text, merge=1)



    if top is None:
        try:
            cmd.save(out_path)
            print("File",out_path + " created\n")
        except:
            print("Error writing file",out_path + '\n')
    else:
        try:
            cmd.save(out_path)
            print("File",out_path + " created\n")
        except:
            print("Error writing file",out_path+'\n')
