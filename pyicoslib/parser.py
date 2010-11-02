"""
Pyicos is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import sys

from lib import argparse
import ConfigParser
from operations import (Turbomix, Extend, Poisson, RemoveRegion, RemoveChromosome, Normalize, Subtract, Trim,
                        Split, Cut, NoWrite, DiscardArtifacts, RemoveDuplicates, OperationFailed, ModFDR, StrandCorrelation)
from core import (BED, ELAND, PK, SPK, ELAND_EXPORT, WIG, CLUSTER_FORMATS, READ_FORMATS, WRITE_FORMATS)
__version__ = '0.8.5'
from defaults import *

class PicosParser:

    def config_section_map(self, section, config_file):
        dict1 = {}
        options = config_file.options(section)
        for option in options:
            try:
                dict1[option] = config_file.get(section, option)
                if dict1[option] == -1:
                    DebugPrint("skip: %s" % option)
            except:
                print("exception on %s!" % option)
                dict1[option] = None
        return dict1


    def new_subparser(self, *args):
        return argparse.ArgumentParser(add_help=False)

    def create_parser(self):
        read_formats = str(READ_FORMATS)
        write_formats = str(WRITE_FORMATS)
        parser = argparse.ArgumentParser(version=__version__)
        subparsers = parser.add_subparsers(help='The operation you want to perform. Note that some operations imply previous automatic operations.')
        #parent parsers
        basic_parser = self.new_subparser()
        basic_parser.add_argument('input', help='The input file or directory. ')
        basic_parser.add_argument('-o','--open-input', action='store_true', default=OPEN_INPUT, help='Defines if the input is half-open or closed notation. [Default %(default)s]')
        basic_parser.add_argument( '-f','--input-format',default=INPUT_FORMAT, help="""The format the input file is written as.
                                 The options are %s. [Default pk]"""%read_formats)
        basic_parser.add_argument('--debug', action='store_true', default=DEBUG)
        basic_parser.add_argument('--no-sort',action='store_true', default=NO_SORT,
                                  help='Force skip the sorting step. WARNING: Use only if you know what you are doing. Processing unsorted files assuming they are will outcome in erroneous results')
        basic_parser.add_argument('--silent' ,action='store_false', default=VERBOSE, dest='verbose', help='Run without printing in screen')
        basic_parser.add_argument('--disable-cache' ,action='store_false', default=CACHED, dest='cached', help='Disable internal reading cache. When Clustering low coverage files, it will increase speed and improve memory usage. With very read dense files, the speed will decrease.')

        output = self.new_subparser()
        output.add_argument('output', help='The output file or directory')

        output_flags = self.new_subparser()
        output_flags.add_argument('-O','--open-output', action='store_true', default=OPEN_OUTPUT, help='Define if the output is half-open or closed notation. [Default %(default)s]')
        output_flags.add_argument('-F','--output-format',default=OUTPUT_FORMAT, help='Format desired for the output. You can choose between %s. WARNING, for some operations, some outputs are not valid. See operation help for more info. [default pk]'%write_formats)

        control = self.new_subparser()
        control.add_argument('control', help='The control file or directory')
        control_format = self.new_subparser()
        control_format.add_argument('--control-format',default=CONTROL_FORMAT, help='The format the control file is written as. [default %(default)s]')
        optional_control = self.new_subparser()
        optional_control.add_argument('--control', help='The control file or directory')
        open_control = self.new_subparser()
        open_control.add_argument('--open-control', action='store_true', default=OPEN_CONTROL, help='Define if the region file is half-open or closed notation. [Default %(default)s]')

        region = self.new_subparser()
        region.add_argument('region', help='The region file')
        optional_region = self.new_subparser()
        optional_region.add_argument('--region', help='The region file or directory')
        region_format = self.new_subparser()
        region_format.add_argument('--region-format',default=REGION_FORMAT, help='The format the region file is written as. [default %(default)s]')
        region_format.add_argument('--open-region', action='store_true', default=OPEN_REGION, help='Define if the region file is half-open or closed notation. [Default %(default)s]')

        label = self.new_subparser()
        label.add_argument('--label', default=LABEL, help='The label that will identify the experiment')

        span = self.new_subparser()
        span.add_argument('--span', default=SPAN, help='The span of the variable and fixed wig formats', type=int)

        round = self.new_subparser()
        round.add_argument('--round',action='store_true',dest='rounding', default=ROUNDING, help='Round the final results to an integer')
        pvalue = self.new_subparser()
        pvalue.add_argument('--p-value',type=float, default=P_VALUE, help='The threshold p-value that will make a cluster significant. [Default %(default)s]')

        tolerated_duplicates =self.new_subparser()
        tolerated_duplicates.add_argument('--duplicates',type=int, default=DUPLICATES, help='The number of duplicated reads accept will be counted. Any duplicated read after this threshold will be discarded. [Default %(default)s]')

        height = self.new_subparser()
        height.add_argument('--height-limit',type=int, default=HEIGHT_LIMIT, help='The cluster height limit Pyicos will analize too. Every cluster that goes up the threshold will have a p-value of 0, therefore considered significant. This parameter is here just for speed purposes, raise it you think that you peak threashold will be over 100 (Almost impossible, but you never know. There is people with crazy data out there.) [Default %(default)s]')

        correction = self.new_subparser()
        correction.add_argument('--correction',type=float, default=CORRECTION, help='This value will correct the size of the genome you are analyzing. This way you can take into consideration the real mappable genome [Default %(default)s]')

        tag_length = self.new_subparser()
        tag_length.add_argument( '--tag-length',default=TAG_LENGTH, type=int, help='The tag length, or the extended one. Needed when converting from a Clustered format (wig, pk) to a non clustered format (bed, eland) [Default %(default)s]')

        frag_size = self.new_subparser()
        frag_size.add_argument('frag_size', help='The estimated inmmunoprecipitated fragment size. This is used by the extend operation to extend the tags, taking into consideration their strand, if provided. If the strand is not provided, Pyicos will assume positive strand.', type=int)
        optional_frag_size = self.new_subparser()
        optional_frag_size.add_argument('-x', '--frag_size', help='The estimated inmmunoprecipitated fragment size. This is used by Pyicos to reconstruct the original signal in the original wet lab experiment.', type=int)

        no_subtract = self.new_subparser()
        no_subtract.add_argument('--no-subtract',action='store_true', default=NO_SUBTRACT, help='Dont subtract the control to the output, only normalize.')

        normalize = self.new_subparser()
        normalize.add_argument('--normalize',action='store_true', default=NORMALIZE, help='Normalize to the control before subtracting')

        split_proportion = self.new_subparser()
        split_proportion.add_argument('--split-proportion', default=SPLIT_PROPORTION, help='Fraction of the cluster height below which the cluster is splitted. [Default %(default)s]', type=float)

        trim_proportion = self.new_subparser()
        trim_proportion.add_argument('--trim-proportion', default=TRIM_PROPORTION, help='Fraction of the cluster height below which the peak is trimmed. Example: For a cluster of height 40, if the flag is 0.05, 40*0.05=2. Every cluster will be trimmed to that height. A position of height 1 is always considered insignificant, no matter what the cluster height is. [Default %(default)s]', type=float)
        trim_absolute = self.new_subparser()
        trim_absolute.add_argument('--trim-absolute', help='The height threshold to trim the clusters.', type=int)

        split_absolute = self.new_subparser()
        split_absolute.add_argument('--split-absolute', default=SPLIT_ABSOLUTE, help='The height threshold to split the clusters. [Default %(default)s]', type=int)

        repeats = self.new_subparser()
        repeats.add_argument('--repeats', help='Number of random repeats when generating the "background" for the modfdr operation[Default %(default)s]', default=REPEATS, type=int)
        masker_file = self.new_subparser()
        masker_file.add_argument('--masker', help='You can provide a masker file that will be used by the modfdr operation background generation so that randomized reads will not fall in this areas')

        discard = self.new_subparser()
        discard.add_argument('--discard', help='Discard the reads that have this particular tag. Example: --discard chr1 will discard all reads with chr1 as tag. You can specify multiple tags to discard using the following notation --discard chr1:chr2:tagN')

        threshold = self.new_subparser()
        threshold.add_argument('--threshold', help='The height threshold used to cut', type=int)

        species = self.new_subparser()
        species.add_argument('-p', '--species', default=SPECIES,
                             help='The species that you are analyzing. This will read the length of the chromosomes of this species from the files inside the folder "chrdesc". If the species information is not known, the filtering step will assume that the chromosomes are as long as the position of the furthest read.')

        correlation_flags = self.new_subparser()
        correlation_flags.add_argument('-x','--max-delta',type=int, default=MAX_DELTA, help='Maximum delta [Default %(default)s]')
        correlation_flags.add_argument('-m','--min-delta',type=int, default=MIN_DELTA, help='Minimum delta [Default %(default)s]')
        correlation_flags.add_argument('--corr-height',type=int, dest='height_filter', default=HEIGHT_FILTER, help='The minimum number of overlapping positive and negative strand reads to include them in the correlation calculation [Default %(default)s]')
        correlation_flags.add_argument('-s','--delta-step',type=int, default=DELTA_STEP, help='The step of the delta values to test [Default %(default)s]')

        protocol_name = self.new_subparser()
        protocol_name.add_argument('protocol_name', help='The protocol configuration file.')
        #callpeaks operation
        subparsers.add_parser('callpeaks', help='The complete peak calling sequence proposed in the future publication. The region file is optional. The same goes for the control file, if not provided, there will not be a normalization or a subtraction.',
                              parents=[basic_parser, optional_control, control_format, open_control, optional_region, output, output_flags, frag_size, round, label, span, no_subtract, discard, pvalue, height, correction, trim_proportion, species, tolerated_duplicates])
        #convert operation
        subparsers.add_parser('convert', help='Convert a file to another file type.',
                              parents=[basic_parser, output, output_flags, round, label, tag_length, span, optional_frag_size])

        #remove chr operation
        parser_chremove = subparsers.add_parser('labelremove', help='Remove all lines that have the specified label(s).', parents=[basic_parser, output, output_flags, round, label])
        parser_chremove.add_argument('discard', help='The tag name (or names) as it appears in the file. Example1: chr1 Example2: chrX:chr3:mytag:myothertag')
        #subtract operation
        subparsers.add_parser('subtract', help='Subtract two pk files. Operating with directories will only give apropiate results if the files and the control are paired in alphabetical order.', parents=[basic_parser, control, control_format, open_control, output, output_flags, round, normalize, tag_length, span, label])
        #split operation
        subparsers.add_parser('split', help='Split the peaks in subpeaks. Only accepts pk or wig as output (other formats under development).', parents=[basic_parser, output, output_flags, round, split_proportion, split_absolute, label])
        #trim operation
        subparsers.add_parser('trim', help='Trim the clusters to a given threshold.', parents=[basic_parser, output, output_flags, round, trim_absolute, label])
        #discard operation
        subparsers.add_parser('discard', help='Discards artifacts from a file. Only accepts pk or wig as output.', parents=[basic_parser, output, output_flags, round, span, label])
        #remove duplicates operation
        subparsers.add_parser('remduplicates', help='Removes the duplicated reads in a file. Only accepts tag files (bed, eland)', parents=[basic_parser, output, output_flags, tolerated_duplicates, round, span, label])
        #normalize operation
        subparsers.add_parser('normalize', help='Normalize a pk file respect of the control.', parents=[basic_parser, control, control_format, output, output_flags, open_control, round, label, span])
        #extend operation
        subparsers.add_parser('extend', help='Extend the reads of a file to the desired length (we currently support only bed and eland files for this operation)', parents=[basic_parser,  output, output_flags, frag_size, round, label, span])
        #poisson analysis
        subparsers.add_parser('poisson', help='Analyze the significance of accumulated reads in the file using the poisson distribution. With this tests you will be able to decide what is the significant threshold for your reads.',
                              parents=[basic_parser, output_flags, frag_size, pvalue, height, correction, species])
        #cut operations
        subparsers.add_parser('filter', help="""Analyze the significance of accumulated reads in the file using the poisson distribution and generate the resulting profiles, in wig or pk formats""",
                              parents=[basic_parser, output, frag_size, output_flags, round, pvalue, height, correction, threshold, species])
        #modfdr analysis
        subparsers.add_parser('modfdr', help="""Use the modified FDR method to determine what clusters are significant in an specific region. Output in a clustered format only.""",
                              parents=[basic_parser, region, output, output_flags, round, pvalue, repeats, masker_file])
        #remove operation
        subparsers.add_parser('remove', help='Removes regions that overlap with another the coordinates in the "black list" file.',
                              parents=[basic_parser, output_flags, region, region_format, output])
        #strcorr operation
        subparsers.add_parser('strcorr', help='A cross-correlation test between forward and reverse strand clusters in order to find the optimal extension length.',
                              parents=[basic_parser, output_flags, correlation_flags])
        #enrichment operation
        #subparsers.add_parser('enrichment', help='An enrichment test', parents=[basic_parser, output_flags, region, region_format, output])
        #protocol reading
        subparsers.add_parser('protocol', help='Import a protocol file to load in Pyicos', parents=[protocol_name])
        return parser

    def run_parser(self):
        parser = self.create_parser()
        parser.set_defaults(input=INPUT, input_format=INPUT_FORMAT, open_input=OPEN_INPUT, debug=DEBUG, discard=DISCARD, output=OUTPUT, control=CONTROL, label = LABEL, output_format=OUTPUT_FORMAT,
                            open_output=OPEN_OUTPUT, rounding=ROUNDING, control_format=CONTROL_FORMAT, region=REGION, region_format=REGION_FORMAT, open_region =OPEN_REGION,
                            frag_size = FRAG_SIZE, tag_length = TAG_LENGTH, span=SPAN, p_value=P_VALUE, height_limit=HEIGHT_LIMIT, correction=CORRECTION, no_subtract = NO_SUBTRACT, normalize = NORMALIZE,
                            trim_proportion=TRIM_PROPORTION,open_control=OPEN_CONTROL, no_sort=NO_SORT, duplicates=DUPLICATES, threshold=THRESHOLD, trim_absolute=TRIM_ABSOLUTE,
                            max_delta=MAX_DELTA, min_delta=MIN_DELTA, height_filter=HEIGHT_FILTER, delta_step=DELTA_STEP, verbose=VERBOSE, species=SPECIES, cached=CACHED, split_proportion=SPLIT_PROPORTION, split_absolute=SPLIT_ABSOLUTE, repeats=REPEATS, masker_file=MASKER_FILE)

        args = parser.parse_args()
        if not args.control_format: #If not specified, the control format is equal to the input format
            args.control_format = args.input_format

        #Add any parameters found in the config file. Override them with anything found in the args later
        if sys.argv[1] == 'protocol':
            config = ConfigParser.ConfigParser()
            config.read(args.protocol_name)
            section = self.config_section_map("Pyicotrocol", config)
            for key, value in section.items(): #this works fine for all string values
                try:
                    t = type(parser._defaults[key])
                    if t == int:
                        args.__dict__[key] = config.getint("Pyicotrocol", key)
                    elif t == float:
                        args.__dict__[key] = config.getfloat("Pyicotrocol", key)
                    elif t == bool:
                        args.__dict__[key] = config.getboolean("Pyicotrocol", key)
                    elif t == str:
                        args.__dict__[key] = config.get("Pyicotrocol", key)

                except KeyError:
                    if key != 'operations':
                        print 'Parameter %s is not a Pyicos parameter'%key

        turbomix = Turbomix(args.input, args.output, args.input_format, args.output_format, args.label, args.open_input, args.open_output, args.debug,
                            args.rounding, args.tag_length, args.discard, args.control, args.control_format, args.open_control, args.region,
                            args.region_format, args.open_region, args.span, args.frag_size, args.p_value, args.height_limit, args.correction,
                            args.trim_proportion, args.no_sort, args.duplicates, args.threshold, args.trim_absolute, args.max_delta,
                            args.min_delta, args.height_filter, args.delta_step, args.verbose, args.species, args.cached, args.split_proportion, args.split_absolute, args.repeats, args.masker_file)


        if sys.argv[1] == 'protocol':
            operations = section['operations'].split(',')
            for operation in operations:
                turbomix.operations.append(operation.strip())
                
        elif sys.argv[1] == 'convert':
            if args.frag_size:
                turbomix.operations = [Extend]

        elif sys.argv[1] == 'subtract':
            turbomix.operations = [Subtract]
            if args.normalize:
                turbomix.operations.append(Normalize)

        elif sys.argv[1] == 'normalize':
            turbomix.operations = [Normalize]

        elif sys.argv[1] == 'extend':
            turbomix.operations = [Extend]

        elif sys.argv[1] == 'strcorr':
            turbomix.operations = [StrandCorrelation, NoWrite]

        elif sys.argv[1] == 'poisson':
            turbomix.operations = [Poisson, NoWrite]

        elif sys.argv[1] == 'filter':
            turbomix.operations = [Poisson, Cut]

        elif sys.argv[1] == 'remove':
            turbomix.operations = [RemoveRegion]

        elif sys.argv[1] == 'chremove':
            turbomix.operations = [RemoveChromosome]

        elif sys.argv[1] == 'split':
            turbomix.operations = [Split]

        elif sys.argv[1] == 'trim':
            turbomix.operations = [Trim]

        elif sys.argv[1] == 'discard':
            turbomix.operations = [DiscardArtifacts]

        elif sys.argv[1] == 'remduplicates':
            turbomix.operations = [RemoveDuplicates]

        elif sys.argv[1] == 'modfdr':
            turbomix.operations = [ModFDR]

        elif sys.argv[1] == 'callpeaks':
            turbomix.operations = [RemoveChromosome, Split, Extend, Poisson, Cut, RemoveDuplicates] 
            if args.duplicates > 1: #If there is only 1 duplicate,  there is no need to discard artifacts
                turbomix.operations.append(DiscardArtifacts)
            if args.region:
                turbomix.operations.append(RemoveRegion)
            if args.control:
                turbomix.operations.append(Normalize)
                if not args.no_subtract:
                    turbomix.operations.append(Subtract)

        
        #parameters are set, now try running
        try:
            turbomix.run()

        except KeyboardInterrupt:
            print 'Canceled by user.'

        except OperationFailed:
            if args.debug:
               raise
            else:
                print 'Operation Failed.'



