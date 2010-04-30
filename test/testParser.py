import unittest
from pyicoslib.parser import PicosParser

class TestParser(unittest.TestCase):
    
    def setUp(self):
        pass

    def test_it_works(self):
        """
        This test checks that the parser will not raise compiling errors.
        """
        try:
            p = PicosParser()
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
