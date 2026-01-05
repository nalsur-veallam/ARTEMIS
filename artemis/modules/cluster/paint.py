from pymol import cmd

it = 0
lastResi = None


colors = [
    "blue",
    "green",
    "red",
    "cyan",
    "yellow",
    "magenta",
    "salmon",
    "lime",
    "hotpink",
    "orange",
    "chartreuse",
    "limegreen",
    "olive",
    "purple",
    "teal",
    "gray",
    "pink",
    "firebrick",
    "chocolate",
    "brown",
    "wheat",
    "violet",
    "aquamarine",
    "palegreen",
    "lightblue",
    "lightpink",
    "skyblue",
    "silver",
    "gold",
]


def draw_clustering(Clusters, out_path, path):
    if Clusters.clustering_labels is None:
        print("Error: clustering labels not found.\n")

    print("\nSCRIPT FOR DISPLAYING CLUSTERS ON THE STRUCTURE IS LAUNCHED\n")

    def coloring(resi, resn, aname, index, b):
        global it, lastResi

        resid = resi + resn

        if not lastResi:
            lastResi = resid

        if lastResi != resid:
            it += 1
            lastResi = resid

        clust = Clusters.clustering_labels[it]

        clust_name = colors[clust - 1] + "_cl_" + str(clust)

        text = (
            "resi "
            + str(resi)
            + " & resn "
            + resn
            + " & name "
            + aname
            + " & index "
            + str(index)
        )
        cmd.select(clust_name, text, merge=1)

        return b

    myspace = {"coloring": coloring}

    cmd.load(path)
    cmd.alter("all", "b = coloring(resi, resn, name, index, b)", space=myspace)
    cmd.show_as("cartoon", cmd.get_object_list("all")[0])

    for i in range(Clusters.NClusters):
        obj = colors[i] + "_cl_" + str(i + 1)
        cmd.color(colors[i], obj)
        cmd.recolor()

    try:
        cmd.save(out_path)
        print("File", out_path + " created\n")
    except Exception as e:
        print("Error writing file", out_path, f" - {e} \n")
