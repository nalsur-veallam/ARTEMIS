import argparse
import numpy as np
from artemis._modules import allostery, convert, _map, cluster

def parse():

      """
      Read command line arguments

      """

      parser = argparse.ArgumentParser()

      subparsers = parser.add_subparsers(title='modules',
                                   description='valid subcommands',
                                   help='ARTEMIS module: map, allostery, cluster (future: convert, analysis)')

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

      map_parser.add_argument('-n0', type=float, default=None,
                              help='Denoised map number of frames (default is inf).')

      map_parser.add_argument('-n1', type=float,
                              help='First custom number of frames.')

      map_parser.add_argument('-n2', type=float,
                              help='Second custom number of frames.')

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
      allostery_parser.add_argument('-verb', action='store_true',
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

      allostery_parser.add_argument('-noseq', type=int, default=0,
                              help='The number of residues that are considered closest in sequence and are not taken into account when calculating the intensity (default is 0).')

      allostery_parser.add_argument('-top', type=float, default=None,
                              help='The number from the interval (0, 100), it filters the top TOP percent when calculating the intensity.')

      allostery_parser.add_argument('-zscore', action='store_true',
                              help='Using z-score format in output.')

      allostery_parser.add_argument('-cluster', type=int, default=None,
                              help='Calculate intensity for custom cluster.')

      allostery_parser.set_defaults(func=allostery)

      # Cluster parser

      cluster_parser = subparsers.add_parser('cluster', help='Clustering module.')

      cluster_parser.add_argument('files', type=str, nargs='+', default='map.json',
                              help=".json-file")

      cluster_parser.add_argument('-verb', action='store_true',
                              help='Write additional output.')

      cluster_parser.add_argument('--cluster', action='store_true',
                              help='Cluster the system according to the mutual information matrix.')

      cluster_parser.add_argument('--draw', action='store_true',
                              help='Draw clusters using PyMol.')

      cluster_parser.add_argument('--analysis', action='store_true',
                              help='Analyze clusters.')

      cluster_parser.add_argument('--study', action='store_true',
                              help='Plot dendrogram and study optimal number of clusters.')

      cluster_parser.add_argument('-nclust', type=int, default=3,
                              help='Number of clusters (default is 3).')

      cluster_parser.add_argument('-o', type=str,
                              help='Outfile path.')

      cluster_parser.add_argument('-strc', type=str,
                              help='Molecule structure path.')

      cluster_parser.add_argument('-min', type=int, default=2,
                              help='Minimum number of clusters for study.')

      cluster_parser.add_argument('-max', type=int, default=10,
                              help='Maximum number of clusters for study.')

      cluster_parser.add_argument('-noseq', type=float, default=0,
                              help='The number of residues that are considered closest in sequence and are not taken into account when calculating the intensity (default is 0).')

      cluster_parser.add_argument('-spectral', action='store_true',
                              help='Use spectral clustering.')

      cluster_parser.set_defaults(func=cluster)

      args = parser.parse_args()

      return args
