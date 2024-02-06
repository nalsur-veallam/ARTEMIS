import argparse
import numpy as np

def parse():

      """
      Read command line arguments

      """

      parser = argparse.ArgumentParser()

      parser.add_argument('--version', action='store_true',
                              help="ARTEMIS version.")


      args = parser.parse_args()

      return args

def parse_map():

      """
      Read command line arguments

      """

      parser = argparse.ArgumentParser()

      parser.add_argument('files', type=str, nargs='+', default='map.json',
                              help=".json-file or .par-file")

      parser.add_argument('--verbose', action='store_true',
                              help='Write additional output.')

      # Draw map runtypes
      parser.add_argument('--draw', action='store_true',
                              help='Draw map')

      parser.add_argument('--denoise', action='store_true',
                              help='Create denoise map')

      parser.add_argument('--gen', action='store_true',
                              help='Generate MI map from .par-file')


      args = parser.parse_args()

      return args

def parse_allostery():

      """
      Read command line arguments

      """

      parser = argparse.ArgumentParser()

      parser.add_argument('files', type=str, nargs='+', default='map.json',
                              help=".json-file")
      parser.add_argument('--verbose', action='store_true',
                              help='Write additional output.')

      parser.add_argument('--search', action='store_true',
                              help='Search allosteric site')

      parser.add_argument('--draw', action='store_true',
                              help='Draw allostery')

      parser.add_argument('--analysis', action='store_true',
                              help='Analyze allostery')


      args = parser.parse_args()

      return args
