from artemis.io.argparser import parse, parse_map, parse_allostery
#from artemis.datatypes import Map
from artemis._modules import allostery, convert, _map
import sys

def main():

    if sys.argv[1] == 'allostery':

        args = parse_allostery()

        allostery(args)

        #raise ValueError('Please select runtype.')

    elif sys.argv[1] == 'map':

        args = parse_map()

        _map(args)

    elif sys.argv[1] == 'convert':

        args = parse_convert(args)

        convert(args)

    else:
        args = parse()



if __name__ == "__main__":
    main()
