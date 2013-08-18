from pyicoteolib.defaults import *
from common import *
import argparse
from ..counter import count_all



def create_parser():
    n_files = new_subparser()
    n_files.add_argument('read_files', nargs='+', help='The files to count from')
    n_files.add_argument("--gtf-file", type=str, help= "The GTF file")
    parser = argparse.ArgumentParser(   version=VERSION, 
                                        description="""Pyicount, part of the Pyicoteo suite. 
                                                     Count read files and generate a pyicos count file for 1 to N read files.""",
                                        parents=[n_files, experiment_flags, output, basic_parser, 
                                                 optional_region, magic_flag],
        )

    return parser

def run_parser(parser, test_args=None):
    args = parser.parse_args()
    print args.read_files

    count_all(args.read_files, args.experiment_format, args.region, args.gtf_file, args.region_magic)


