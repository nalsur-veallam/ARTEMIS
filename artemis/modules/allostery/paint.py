import numpy as np
import pandas as pd
from pymol import cmd

k = 0


def max_top(array, top, NResidues):
    size = NResidues
    top10 = int(top / 100 * size)

    supp = np.zeros(size)
    for i in range(top10):
        supp[i] = 1

    d = {"data": array, "index": np.arange(0, NResidues)}

    df = pd.DataFrame(data=d)
    df = df.sort_values(by=["data"], ascending=False)

    for i in range(NResidues):
        if np.array(df["data"])[i] < 0:
            supp[i] = 1

    df["data"] = df["data"] * supp
    df = df.sort_values(by=["index"])
    return np.array(df["data"])


def draw_allostery(Allostery, out_path, top, noseq, str_path):
    visual = "Y"
    act = True
    alls = True

    print(
        "\nSCRIPT FOR DISPLAYING ALLOSTERIC COMMUNICATION INTENSITY ON THE STRUCTURE IS LAUNCHED\n"
    )

    if Allostery.allosteric_site is None:
        alls = False
    else:
        alls_arr = []
        group_alls = []
        group_alls_names = []
        group_alls_chains = []

        for i in range(Allostery.NResidues):
            if i + 1 in Allostery.allosteric_site:
                alls_arr.append(1)
            else:
                alls_arr.append(0)

    intensity = np.array(Allostery.intensity)
    if top is not None:
        intensity = max_top(intensity, top, Allostery.NResidues)

    for i in range(Allostery.NResidues):
        if i + 1 in Allostery.active_site:
            intensity[i] = -1
    intensity[intensity < 0] = -1

    resids = []
    group = []
    group_names = []
    group_chains = []

    def renumbering(chain, resi, resn, b):
        resid = chain + resi + resn
        if resid not in resids:
            resids.append(resid)
            global k
            k += 1

            if intensity[k - 1] < 0:
                group.append(resi)
                group_names.append(resn)
                group_chains.append(chain)

            if alls:
                if alls_arr[k - 1] == 1:
                    group_alls.append(resi)
                    group_alls_names.append(resn)
                    group_alls_chains.append(chain)

            return intensity[k - 1]
        return intensity[k - 1]

    myspace = {"renumbering": renumbering}

    cmd.load(str_path)
    cmd.alter("all", "b = renumbering(chain, resi, resn, b)", space=myspace)

    if visual == "Y":
        obj = cmd.get_object_list("all")[0]
        cmd.show_as("cartoon", obj)
        cmd.cartoon("putty", obj)
        cmd.set("cartoon_putty_scale_min", min(intensity), obj)
        cmd.set("cartoon_putty_scale_max", max(intensity), obj)
        cmd.set("cartoon_putty_transform", 0, obj)
        cmd.set("cartoon_putty_radius", 0.4, obj)
        cmd.spectrum("b", "rainbow", obj)
        cmd.ramp_new("count", obj, [min(intensity), max(intensity)], "rainbow")
        cmd.recolor()

    for j in range(len(group)):
        text = (
            "resi "
            + str(group[j])
            + "& resn "
            + str(group_names[j])
            + "& chain "
            + str(group_chains[j])
        )
        cmd.select("active_site", text, merge=1)

    if act:
        obj = "active_site"
        cmd.show_as("spheres", obj)
        cmd.color("white", obj)
        cmd.recolor()

    if alls:
        for j in range(len(group_alls)):
            text = (
                "resi "
                + str(group_alls[j])
                + "& resn "
                + str(group_alls_names[j])
                + "& chain "
                + str(group_alls_chains[j])
            )
            cmd.select("allosteric_site", text, merge=1)

    try:
        cmd.save(out_path)
        print("File", out_path + " created\n")
    except Exception as e:
        print(f"Error writing file {out_path}: {e}\n")
