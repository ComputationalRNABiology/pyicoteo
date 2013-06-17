import sys
import unittest

from pyicoteolib.parser.utils import set_defaults, init_turbomix, parse_validate_args
from pyicoteolib.parser import pyicoclip, pyicoenrich, pyicoller, pyicoregion, pyicos


'''
Sometimes it may be useful to have an ArgumentParser parse args other than those of sys.argv. 
This can be accomplished by passing a list of strings to parse_args. This is useful for testing at the interactive prompt:



>>> parser = argparse.ArgumentParser()
>>> parser.add_argument(
...     'integers', metavar='int', type=int, choices=xrange(10),
...  nargs='+', help='an integer in the range 0..9')
>>> parser.add_argument(
...     '--sum', dest='accumulate', action='store_const', const=sum,
...   default=max, help='sum the integers (default: find the max)')



>>> parser.parse_args(['1', '2', '3', '4'])
Namespace(accumulate=<built-in function max>, integers=[1, 2, 3, 4])
>>> parser.parse_args('1 2 3 4 --sum'.split())
Namespace(accumulate=<built-in function sum>, integers=[1, 2, 3, 4])
'''

TEST_BED1 = "test_files/mini_sorted.bed"
TEST_BED2 = "test_files/mini_sorted2.bed"
TEST_SAM = "test_files/p300_sample.sam"
REGION = "test_files/region_test.bed"
GTF = "test_files/head500_gencode.gtf"
RESULTS_DIR = "test_files/results"


class TestParser(unittest.TestCase):
     

    def test_pyicoclip(self):
        parser = pyicoclip.create_parser()
        pyicoclip.run_parser(parser, ("%s %s/clipregion_out --region %s -f sam --silent"%(TEST_SAM, RESULTS_DIR, REGION)).split())
        pyicoclip.run_parser(parser, ("%s %s/clipmagic_out -f sam --gff-file %s --region-magic genebody --silent"%(TEST_SAM, REGION, RESULTS_DIR)).split())


    def test_pyicoregion(self):
        parser = pyicoregion.create_parser()
        pyicoregion.run_parser(parser, ("%s %s/pyicoregion_out --region-magic exons --silent"%(GTF, RESULTS_DIR)).split())

    def test_pyicoller(self):
        parser = pyicoller.create_parser()
        pyicoller.run_parser(parser, ("%s %s/pyicoller_out -f bed --silent"%(TEST_BED1, RESULTS_DIR)).split())

    def test_pyicoenrich(self):
        parser = pyicoenrich.create_parser()
        pyicoenrich.run_parser(parser, ("-reads %s %s -output %s/enrich_out -f bed --silent"%(TEST_BED1, TEST_BED2, RESULTS_DIR)).split())

    def test_pyicos(self):
        #TODO iterate sub-commands
        parser = pyicos.create_parser()       
        sys.argv.append("subtract") #FIXME workaround, ideally pyicos parser should not use the sys library
        pyicos.run_parser(parser, ("subtract %s %s %s/subtract_out -f bed --silent"%(TEST_BED1, TEST_BED2, RESULTS_DIR)).split())

def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestParser))
    return suite


if __name__ == '__main__':
    unittest.main()
