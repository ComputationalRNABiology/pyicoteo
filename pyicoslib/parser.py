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
        basic_parser.add_argument('-o','--open-input', action='store_true', default=False, help='Defines if the input is half-open or closed notation. [Default %(default)s]')
        basic_parser.add_argument( '-f','--input-format',default='pk', help="""The format the input file is written as.
                                 The options are %s. [Default pk]"""%read_formats)
        basic_parser.add_argument('--debug', action='store_true', default=False)
        basic_parser.add_argument('--no-sort',action='store_true', default=False,
                                  help='Force skip the sorting step. WARNING: Use only if you know what you are doing. Processing unsorted files assuming they are will outcome in erroneous results')
        basic_parser.add_argument('--silent' ,action='store_false', default=True, dest='verbose', help='Force skip the sorting step. WARNING: Working with unsorted files will outcome in unexpected results')

        output = self.new_subparser()
        output.add_argument('output', help='The output file or directory')

        output_flags = self.new_subparser()
        output_flags.add_argument('-O','--open-output', action='store_true', default=False, help='Define if the output is half-open or closed notation. [Default %(default)s]')
        output_flags.add_argument('-F','--output-format',default='pk', help='Format desired for the output. You can choose between %s. WARNING, for some operations, some outputs are not valid. See operation help for more info. [default pk]'%write_formats)

        control = self.new_subparser()
        control.add_argument('control', help='The control file or directory')
        control_format = self.new_subparser()
        control_format.add_argument('--control-format',default='pk', help='The format the control file is written as. [default %(default)s]')
        optional_control = self.new_subparser()
        optional_control.add_argument('--control', help='The control file or directory')
        open_control = self.new_subparser()
        open_control.add_argument('--open-control', action='store_true', default=False, help='Define if the region file is half-open or closed notation. [Default %(default)s]')

        region = self.new_subparser()
        region.add_argument('region', help='The region file')
        optional_region = self.new_subparser()
        optional_region.add_argument('--region', help='The region file or directory')
        region_format = self.new_subparser()
        region_format.add_argument('--region-format',default='bed', help='The format the region file is written as. [default %(default)s]')
        region_format.add_argument('--open-region', action='store_true', default=False, help='Define if the region file is half-open or closed notation. [Default %(default)s]')

        label = self.new_subparser()
        label.add_argument('--label', default='pyicos_output', help='The label that will identify the experiment')

        span = self.new_subparser()
        span.add_argument('--span', default=25, help='The span of the variable and fixed wig formats', type=int)

        round = self.new_subparser()
        round.add_argument('--round',action='store_true',dest='rounding', default=False, help='Round the final results to an integer')
        pvalue = self.new_subparser()
        pvalue.add_argument('--p-value',type=float, default=0.01, help='The p-value we consider to be significant in our statistical test. [Default %(default)s]')

        tolerated_duplicates =self.new_subparser()
        tolerated_duplicates.add_argument('--duplicates',type=int, default=3, help='The number of duplicates we accept as valid. ()[Default %(default)s]')

        height = self.new_subparser()
        height.add_argument('--height-limit',type=int, default=100,help='The cluster height limit Pyicos will analize too. Every cluster that goes up the threshold will have a p-value of 0, therefore considered significant. This parameter is here just for speed purposes, raise it you think that you peak threashold will be over 100 (Almost impossible, but you never know. There is people with crazy data out there.) [Default %(default)s]')

        correction = self.new_subparser()
        correction.add_argument('--correction',type=float, default=1., help='This value will correct the size of the genome you are analyzing. This way you can take into consideration the real mappable genome [Default %(default)s]')

        tag_length = self.new_subparser()
        tag_length.add_argument( '--tag-length',default=None, type=int, help='The tag length, or the extended one. Needed when converting from a Clustered format (wig, pk) to a non clustered format (bed, eland) [Default %(default)s]')

        frag_size = self.new_subparser()
        frag_size.add_argument('frag_size', help='The estimated inmmunoprecipitated fragment size. This is used by Pyicos to reconstruct the original signal in the original wet lab experiment.', type=int)
        optional_frag_size = self.new_subparser()
        optional_frag_size.add_argument('-x', '--frag_size', dest='frag_size', help='The estimated inmmunoprecipitated fragment size. This flag is optional.', type=int)

        no_subtract = self.new_subparser()
        no_subtract.add_argument('--no-subtract',action='store_true', default=False, help='Dont subtract the control to the output, only normalize.')

        normalize = self.new_subparser()
        normalize.add_argument('--normalize',action='store_true', default=False, help='Normalize to the control before subtracting')

        trim_proportion = self.new_subparser()
        trim_proportion.add_argument('--trim-proportion', default=0.3, help='Fraction of the cluster height below which the peak is trimmed. Example: For a cluster of height 40, if the flag is 0.05, 40*0.05=2. Every cluster will be trimmed to that height. A position of height 1 is always considered insignificant, no matter what the cluster height is. [Default %(default)s]', type=float)
        trim_absolute = self.new_subparser()
        trim_absolute.add_argument('--trim-absolute', help='The height threshold used to split or trim the clusters.', type=int)

        discard = self.new_subparser()
        discard.add_argument('--discard', help='Discard the reads that have this particular tag. Example: --discard chr1 will discard all reads with chr1 as tag. You can specify multiple tags to discard using the following notation --discard chr1:chr2:tagN')

        threshold = self.new_subparser()
        threshold.add_argument('--threshold', help='The height threshold used to cut', type=int)

        correlation_flags = self.new_subparser()
        correlation_flags.add_argument('-x','--max-delta',type=int, default=500, help='Maximum delta [Default %(default)s]')
        correlation_flags.add_argument('-m','--min-delta',type=int, default=0, help='Minimum delta [Default %(default)s]')
        correlation_flags.add_argument('-t','--height-filter',type=int, default=8, help='Height to filter the peaks [Default %(default)s]')
        correlation_flags.add_argument('-s','--delta-step',type=int, default=1, help='The step of the delta values to test [Default %(default)s]')

        protocol_name = self.new_subparser()
        protocol_name.add_argument('protocol_name', help='The protocol configuration file.')
        #callpeaks operation
        subparsers.add_parser('callpeaks', help='The complete peak calling sequence proposed in the future publication. The region file is optional. The same goes for the control file, if not provided, there will not be a normalization or a subtraction.', parents=[basic_parser, optional_control, control_format, open_control, optional_region, output, output_flags, frag_size, round, label, span, no_subtract, discard, pvalue, height, correction, trim_proportion])
        #convert operation
        subparsers.add_parser('convert', help='Convert a file to another file type.', parents=[basic_parser,  output, output_flags, round, label, tag_length, span, optional_frag_size])
        #remove chr operation
        parser_chremove = subparsers.add_parser('labelremove', help='Remove all lines that have the specified label(s).', parents=[basic_parser, output, output_flags, round, label])
        parser_chremove.add_argument('discard', help='The tag name (or names) as it appears in the file. Example1: chr1 Example2: chrX:chr3:mytag:myothertag')
        #subtract operation
        subparsers.add_parser('subtract', help='Subtract two pk files. Operating with directories will only give apropiate results if the files and the control are paired in alphabetical order.', parents=[basic_parser, control, control_format, open_control, output, output_flags, round, normalize, tag_length, span, label])
        #split operation
        subparsers.add_parser('split', help='Split the peaks in subpeaks. Only accepts pk or wig as output (other formats under development).', parents=[basic_parser, output, output_flags, round, trim_proportion, trim_absolute, label])
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
                              parents=[basic_parser, output_flags, frag_size, pvalue, height, correction])
        #cut operations
        subparsers.add_parser('filter', help="""Analyze the significance of accumulated reads in the file using the poisson distribution and generate the resulting profiles, in wig or pk formats""",
                              parents=[basic_parser, output, frag_size, output_flags, round, pvalue, height, correction, threshold])
        #modfdr analysis
        subparsers.add_parser('modfdr', help="""Use the modified FDR method to determine what clusters are significant in an specific region. Output in a clustered format only.""",
                              parents=[basic_parser, region, output, output_flags, round])
        #remove operation
        subparsers.add_parser('remove', help='Removes regions that overlap with another the coordinates in the "black list" file.',
                              parents=[basic_parser, output_flags, region, region_format, output])
        #strcorr operation
        subparsers.add_parser('strcorr', help='A cross-correlation test between forward and reverse strand clusters in order to find the optimal extension length.',
                              parents=[basic_parser, output_flags, correlation_flags])
        #strcorr operation
        subparsers.add_parser('protocol', help='Import a protocol file to load in Pyicos', parents=[protocol_name])
        return parser

    def __init__(self):
        parser = self.create_parser()
        parser.set_defaults(input='', input_format=PK, open_input=False, debug=False, discard=0, output='', control='', label = 'noname', output_format=PK,
                            open_output =False, rounding=False, control_format=PK, region='', region_format=BED, open_region = False,
                            frag_size = 0, tag_length = 0, span=40, p_value=0.01, height_limit=100, correction=1., no_subtract = False, normalize = False,
                            trim_proportion=0.3,open_control=False, no_sort=False, duplicates=3, threshold=0, trim_absolute=0,
                            max_delta=500, min_delta=0, height_filter=8, delta_step=1, verbose=True)


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
                            args.min_delta, args.height_filter, args.delta_step, args.verbose)


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
            turbomix.operations = [RemoveChromosome, Split, Extend, DiscardArtifacts, Poisson, Cut]
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



