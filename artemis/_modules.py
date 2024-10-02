from artemis.modules.map.draw import draw_map
from artemis.modules.map.contacts import contacts_map
from artemis.datatypes import MAP, GROUPS, ALLOSTERY, CLUSTERS
import artemis
import os
import subprocess
from artemis.modules.allostery.search import search_allostery
from artemis.modules.allostery.paint import draw_allostery
from artemis.modules.allostery.analysis import analyze_allostery
from artemis.modules.cluster.clustering import do_clustering
from artemis.modules.cluster.clustering import do_spectral_clustering
from artemis.modules.cluster.paint import draw_clustering
from artemis.modules.cluster.study import study_clustering
from artemis.modules.cluster.analysis import analyze_clustering
from artemis.modules.analysis import mimsa
def allostery(args):

    for f in args.files:
        Allostery = ALLOSTERY(f)
        if Allostery.is_ok():
            break
    if not Allostery.is_ok():

        for f in args.files:
            Map = MAP(f)
            if Map.is_ok():
                break
        if not Map.is_ok():
            Allostery.interrupt()
        else:
            for f in args.files:
                Groups = GROUPS(f)
                if Groups.is_ok() and "active_site" in Groups.data.keys():
                    break
            if Groups.is_ok() and "active_site" in Groups.data.keys():
                Allostery = ALLOSTERY()
                Allostery.from_map_and_groups(Map, Groups)

                if not Allostery.is_ok():
                    Allostery.interrupt()
            else:
                Allostery.interrupt()

    if args.search:

        if args.o is None:
            args.o = 'allosteric_intensity.pdf'

        if args.cluster is not None:
            calculated = False

            for f in args.files:
                Clusters = CLUSTERS(f)
                if Clusters.is_ok() and Clusters.exist() and Clusters.clustering_labels is not None:
                    search_allostery(Allostery, args.o, args.top, args.noseq, args.zscore, Clusters.clustering_labels, args.cluster)
                    calculated = True
                    break
                else:
                    pass

            if not calculated:
                search_allostery(Allostery, args.o, args.top, args.noseq, args.zscore)

        else:
            search_allostery(Allostery, args.o, args.top, args.noseq, args.zscore)

    if args.draw:

        if args.o is None:
            args.o = 'allosteric_intensity.pse'

        draw_allostery(Allostery, args.o, args.top, args.noseq, args.strc)

    if args.analysis:
        if args.o is None:
            args.o = 'allosteric_analysis.pdf'

        analyze_allostery(Allostery, args.o, args.top, args.noseq, args.zscore)


def _map(args):

    if args.denoise:
        path = os.path.dirname(artemis.__file__)

        if args.o is None:
            args.o = 'map.json'

        if args.lin:
            if args.n0 is None:
                subprocess.run(path+"/../bin/denoise -o " + str(args.o) + " -lin -n1 " + str(args.n1) + " -n2 " + str(args.n2) + " -f1 " + args.files[0] + " -f2 " + args.files[1], shell=True)
            else:
                subprocess.run(path+"/../bin/denoise -o " + str(args.o) + " -lin -n0 " + str(args.n0) + " -n1 " + str(args.n1) + " -n2 " + str(args.n2) + " -f1 " + args.files[0] + " -f2 " + args.files[1], shell=True)
        else:
            if args.n0 is None:
                subprocess.run(path+"/../bin/denoise -o " + str(args.o) + " -n1 " + str(args.n1) + " -n2 " + str(args.n2) + " -f1 " + args.files[0] + " -f2 " + args.files[1], shell=True)
            else:
                subprocess.run(path+"/../bin/denoise -o " + str(args.o) + " -n0 " + str(args.n0) + " -n1 " + str(args.n1) + " -n2 " + str(args.n2) + " -f1 " + args.files[0] + " -f2 " + args.files[1], shell=True)

    elif args.gen:
        path = os.path.dirname(artemis.__file__)

        if args.o is None:
            args.o = 'map.json'

        subprocess.run(path+"/../bin/get_map -o " + str(args.o) + " -f " + str(args.files[0]), shell=True)

    elif args.draw:

        if args.o is None:
            args.o = 'map.pdf'

        for f in args.files:
            Map = MAP(f)
            if Map.is_ok():
                break
        if not Map.is_ok():
            Map.interrupt()

        draw_map(Map, args.o, args.diag, args.norm)

    elif args.contacts:

        if args.o is None:
            args.o = 'contacs.pdf'

        for f in args.files:
            Map = MAP(f)
            if Map.is_ok():
                break
        if not Map.is_ok():
            Map.interrupt()

        contacts_map(Map, args.o, args.diag, args.norm, args.vmax)

def cluster(args):
    for f in args.files:
        Clusters = CLUSTERS(f, NClusters=args.nclust)
        if Clusters.is_ok() and Clusters.exist():
            for f in args.files:
                Groups = GROUPS(f)
                if Groups.is_ok():
                     Clusters.read_groups(Groups)
            break

    if not (Clusters.is_ok() and Clusters.exist()):

        for f in args.files:
            Map = MAP(f)
            if Map.is_ok():
                break
        if not Map.is_ok():
            Clusters.interrupt()
        else:
            for f in args.files:
                Groups = GROUPS(f)
                if Groups.is_ok() and "reference_group" in Groups.data.keys():
                    break
            if Groups.is_ok() and "restriction_group" in Groups.data.keys():
                Clusters = CLUSTERS()
                Clusters.from_map_and_groups(Map, Groups)
            else:
                Clusters = CLUSTERS(NClusters=args.nclust)
                Clusters.from_map_and_groups(Map, Groups)

            if not Clusters.is_ok():
                Clusters.interrupt()

    if args.cluster:

        if args.spectral:
            if args.o is None:
                args.o = 'spectral_clustering.pdf'

            do_spectral_clustering(Clusters, args.o)

        else:

            if args.o is None:
                args.o = 'clustering.pdf'

            do_clustering(Clusters, args.o)

    if args.draw:

        if args.o is None:
            args.o = 'clustering.pse'

        draw_clustering(Clusters, args.o, args.strc)

    if args.analysis:

        if args.o is None:
            args.o = 'clusters_analysis.pdf'

        analyze_clustering(Clusters, args.noseq, args.o)

    if args.study:

        if args.o is None:
            args.o = 'clustering.pdf'

        study_clustering(Clusters, args.min, args.max, args.o, args.spectral)

def analysis(args):
    if args.mimsa:
        mimsa.map_from_msa(args.n, args.o, args.f, args.apc, args.rcw, args.cl, args.igg, args.zs, args.igc,
                           args.igp, args.ign)
def convert(args):
    pass
