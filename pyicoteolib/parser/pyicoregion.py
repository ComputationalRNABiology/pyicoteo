from pyicoteolib.defaults import *
from utils import *
import argparse

def create_parser():
    region_parser = new_subparser()
    region_parser.add_argument('--gff-file', help="GFF input file")
    output = new_subparser()
    output.add_argument('-output', help='The output file')
    region_parser.add_argument('--region-magic', nargs='+', help="Desired features to filter (exons, introns, sliding window for inter-/intragenic zones, tss)")



    parser = argparse.ArgumentParser(version=VERSION, 
                                     description='Standalone region operations', 
                                     parents=[region_parser, output]
                                     )
    return parser

def run_parser():
    parser = create_parser()
    args = parse_validate_args(parser)

    turbomix = init_turbomix(args, parser_name="pyicoregion")

    turbomix.operations = [REGIONS]

    turbomix.run()


