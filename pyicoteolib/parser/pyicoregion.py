from pyicoteolib.defaults import *
from utils import *
import argparse

def create_parser():
    region_parser = new_subparser()
    region_parser.add_argument('gff-file', help="GFF input file")
    output = new_subparser()
    output.add_argument('output', help='The output file')

    parser = argparse.ArgumentParser(version=VERSION, 
                                     description='Standalone region operations', 
                                     parents=[region_parser, output, output_flags, basic_parser, magic_flag]
                                     )
    return parser

def run_parser(parser, test_args=None):
    args = parse_validate_args(parser, test_args)
    turbomix = init_turbomix(args, parser_name="pyicoregion")
    turbomix.operations = [REGIONS]
    turbomix.run()


