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

class TestParser(unittest.TestCase):
    
    def setUp(self):
        pass

    def parser_test(self, my_module):
        """
        Tests that the parsers can be created 
        """
        try:
            parser = my_module.create_parser()

        except SystemExit, e:
            self.assertEquals(type(e), type(SystemExit()))
            self.assertEquals(e.code, 2)

        except Exception, e:
            self.fail("Unexpected Exception %s"%e)

    def test_pyicoclip(self):
        self.parser_test(pyicoclip)

    def test_pyicoregion(self):
        self.parser_test(pyicoregion)

    def test_pyicoller(self):
        self.parser_test(pyicoller)

    def test_pyicoenrich(self):
        #self.parser_test(pyicoenrich)
        parser = pyicoenrich.create_parser()
        pyicoenrich.run_parser(parser, "-reads test_files/mini_sorted.bed test_files/mini_sorted2.bed -output test_files/results/enrich_out -f bed --silent".split())

        #args = parser.parse_args('-reads a_sample.bed b_sample.bed -output bla.txt -f bed'.split())



    def test_pyicos(self):
        #TODO iterate sub-commands
        self.parser_test(pyicos)        


def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestParser))
    return suite


if __name__ == '__main__':
    unittest.main()
