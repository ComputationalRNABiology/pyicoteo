import unittest
from pyicoslib.utils import DualSortedReader
from pyicoslib.core import BED

class TestUtils(unittest.TestCase):
    
    def test_dual_reader(self):
        reader = DualSortedReader("test_files/mini_sorted.bed", "test_files/mini_sorted2.bed", BED)
        for line in reader:
            print line,

def suite():
   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(TestUtils))
   return suite
