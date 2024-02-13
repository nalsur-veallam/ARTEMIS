from artemis.modules.map.draw import draw_map
from artemis.datatypes import MAP, GROUPS, ALLOSTERY
import artemis
import os
import subprocess
from artemis.modules.allostery.search import search_allostery
from artemis.modules.allostery.paint import draw_allostery

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

        search_allostery(Allostery, args.o, args.top, args.noseq, args.zscore)

    if args.draw:

        if args.o is None:
            args.o = 'allosteric_intensity.pse'

        draw_allostery(Allostery, args.o, args.top, args.noseq, args.strc)



    pass

def _map(args):
    if args.draw:

        if args.o is None:
            args.o = 'map.pdf'

        for f in args.files:
            Map = MAP(f)
            if Map.is_ok():
                break
        if not Map.is_ok():
            Map.interrupt()

        draw_map(Map, args.o, args.diag, args.norm)
    elif args.denoise:
        path = os.path.dirname(artemis.__file__)

        if args.o is None:
            args.o = 'map.json'

        if args.lin:
            subprocess.run(path+"/../bin/denoise -o " + str(args.o) + " -lin -dt0 " + str(args.dt0) + " -dt1 " + str(args.dt1) + " -dt2 " + str(args.dt2) + " -f1 " + args.files[0] + " -f2 " + args.files[1], shell=True)
        else:
            subprocess.run(path+"/../bin/denoise -o " + str(args.o) + " -dt0 " + str(args.dt0) + " -dt1 " + str(args.dt1) + " -dt2 " + str(args.dt2) + " -f1 " + args.files[0] + " -f2 " + args.files[1], shell=True)

    elif args.gen:
        path = os.path.dirname(artemis.__file__)

        if args.o is None:
            args.o = 'map.json'

        subprocess.run(path+"/../bin/get_map -o " + str(args.o) + " -f " + str(args.files[0]), shell=True)


def convert(args):
    pass
