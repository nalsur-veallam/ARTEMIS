from artemis.io.argparser import parse, parse_map, parse_allostery
from artemis.datatypes import Map
from artemis.map.map import map
from artemis.allostery.draw import draw
import sys

def main():

    if sys.argv[1] == 'allostery':
        args = parse_allostery()
        if args.search:
            search_allostery(args)
        else:
            raise ValueError('Please select runtype.')

    elif sys.argv[1] == 'map':
        args = parse_map()
        if args.draw:
            draw_map(args)
        else:
            raise ValueError('Please select runtype.')
    else:
        args = parse()


if __name__ == "__main__":
    main()
