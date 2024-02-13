import argparse
import numpy as np
from artemis._modules import allostery, convert, _map

def parse():

      """
      Read command line arguments

      """

      parser = argparse.ArgumentParser()

      subparsers = parser.add_subparsers(title='modules',
                                   description='valid subcommands',
                                   help='ARTEMIS module: map, allostery (future: cluster, convert, analysis)')

      parser.add_argument('--version', action='store_true',
                              help="ARTEMIS version.")

      #Map parser

      map_parser = subparsers.add_parser('map', help='MIE map module.')

      map_parser.add_argument('files', type=str, nargs='+', default='map.json',
                              help=".json-file or .par-file")

      map_parser.add_argument('-verb', action='store_true',
                              help='Write additional output.')

      map_parser.add_argument('--draw', action='store_true',
                              help='Draw map')

      map_parser.add_argument('--denoise', action='store_true',
                              help='Create denoised map.')

      map_parser.add_argument('--gen', action='store_true',
                              help='Generate MI map (json) from .par-file.')

      map_parser.add_argument('-dt0', type=float, default=0,
                              help='Denoised map dt (default is 0).')

      map_parser.add_argument('-dt1', type=float,
                              help='First custom map dt.')

      map_parser.add_argument('-dt2', type=float,
                              help='Second custom map dt.')

      map_parser.add_argument('-o', type=str,
                              help='Outfile path.')

      map_parser.add_argument('-lin', action='store_true',
                              help='Linear denoising.')

      map_parser.add_argument('-norm', action='store_true',
                              help='Normalize the MIE matrix before drawing.')

      map_parser.add_argument('-diag', action='store_true',
                              help='Draw the diagonal of the MIE matrix.')

      map_parser.set_defaults(func=_map)

      # Allostery parser

      allostery_parser = subparsers.add_parser('allostery', help='Allostery module.')

      allostery_parser.add_argument('files', type=str, nargs='+', default='map.json',
                              help=".json-file")
      allostery_parser.add_argument('--verb', action='store_true',
                              help='Write additional output.')

      allostery_parser.add_argument('--search', action='store_true',
                              help='Search allosteric site.')

      allostery_parser.add_argument('--draw', action='store_true',
                              help='Draw allostery using PyMol.')

      allostery_parser.add_argument('--analysis', action='store_true',
                              help='Analyze allostery.')

      allostery_parser.add_argument('-o', type=str,
                              help='Outfile path.')

      allostery_parser.add_argument('-strc', type=str,
                              help='Molecule structure path.')

      allostery_parser.add_argument('-noseq', type=float, default=0,
                              help='The number of residues that are considered closest in sequence and are not taken into account when calculating the intensity (default is 0).')

      allostery_parser.add_argument('-top', type=float,
                              help='The number from the interval (0, 100), it filters the top TOP percent when calculating the intensity.')

      allostery_parser.add_argument('-zscore', action='store_true',
                              help='Using z-score format in output.')

      allostery_parser.set_defaults(func=allostery)

      args = parser.parse_args()

      return args
