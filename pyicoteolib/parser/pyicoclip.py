from pyicoteolib.defaults import *
from utils import *
import argparse

def create_parser():
    parser = argparse.ArgumentParser(version=VERSION, 
                                     description="Pyicoclip, part of the Pyicoteo suite, is a peak caller specifically designed for CLIP-Seq analysis, based on the ModFDR method proposed by Yeo et al.",
                                     parents=[experiment, experiment_flags, basic_parser, region, output, output_flags, round, pvalue, repeats, remlabels]
                                     )
    return parser

def run_parser(parser):
    args = parse_validate_args(parser)
    turbomix = init_turbomix(args, parser_name="pyicoclip")
    turbomix.operations = [MODFDR]
    turbomix.run()


