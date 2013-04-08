import unittest
from pyicoteolib.parser import PicosParser

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

    def test_it_works(self):
        """
        This test checks that the parser will not raise compiling errors.
        """
        try:
            p = PicosParser()
            p.run_parser()
        except SystemExit, e:
            self.assertEquals(type(e), type(SystemExit()))
            self.assertEquals(e.code, 2)
        except Exception, e:
            self.fail('unexpected exception: %s' % e)
        else:
            self.fail('SystemExit exception expected')

        

def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestParser))
    return suite


if __name__ == '__main__':
    unittest.main()
