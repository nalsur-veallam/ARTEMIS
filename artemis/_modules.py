from artemis.modules.map.draw import draw_map
import subprocess
#from artemis.modules.allostery.search import draw_allostery

def allostery(args):
    pass

def _map(args):
    if args.draw:
        draw_map(args)
    elif args.denoise:
        if args.lin:
            path = subprocess.run("pwd artemis", shell=True, text=True, capture_output=True)
            subprocess.run(path.stdout[:-1]+"/bin/denoise -o map.json -lin -n map -dt1 " + str(args.dt1) + " -dt2 " + str(args.dt2) + " -f1 " + args.files[1] + " -f2 " + args.files[2], shell=True)


def convert(args):
    pass
