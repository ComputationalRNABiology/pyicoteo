from pyicoteolib.defaults import *
from utils import *
import argparse

def create_parser():
    parser = argparse.ArgumentParser(version=VERSION, 
                                     description='An enrichment test based on the MA plots using mapped reads files. Pyicoenrich output will consist in a results table and a MA plot (optional, but matplotlib required >=0.9.7). The fields of this table are as follows: %s'%(" | ".join(enrichment_keys)), 
                                     parents=[exp_or_count, experiment_flags, basic_parser, output_flags, optional_replica, optional_region, 
                                              region_format, optional_output, enrichment_flags, tmm_flag, quant_flag, total_reads_flags, 
                                              pseudocount, zscore]
                                     )
    return parser

def run_parser():
    parser = create_parser()
    args = parse_validate_args(parser)

    if args.counts_file: #the formats are overridden when using enrichment (only of cosmetic value, when printing the flags)   
        args.experiment_format = COUNTS
        args.experiment_b_format = COUNTS
        args.output_format = COUNTS

    if not args.control_format: #If not specified, the control format is equal to the experiment format
        args.control_format = args.experiment_format
        args.open_control = args.open_experiment

    if args.experiments:
        args.experiment, args.experiment_b = args.experiments

    turbomix = init_turbomix(args, parser_name="pyicoenrich")


    turbomix.operations = [ENRICHMENT, CALCZSCORE]
    if not args.skip_plot:
       turbomix.operations.append(PLOT)
    turbomix.run()


