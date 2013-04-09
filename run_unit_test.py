#!/usr/bin/env python2

import unittest

from utest import testCore, testOperation, testUtils

suite1 = testCore.suite()
suite2 = testOperation.suite()
#suite3 = testParser.suite()
suite4 = testUtils.suite()

suite = unittest.TestSuite()
suite.addTest(suite1)
suite.addTest(suite2)
#suite.addTest(suite3)
suite.addTest(suite4)

unittest.TextTestRunner(verbosity=2).run(suite)


