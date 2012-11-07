#!/usr/bin/env python

import unittest

from test import testCore, testOperation, testParser, testUtils

suite1 = testCore.suite()
suite2 = testOperation.suite()
suite3 = testParser.suite()
suite4 = testUtils.suite()

suite = unittest.TestSuite()
suite.addTest(suite1)
suite.addTest(suite2)
suite.addTest(suite3)
suite.addTest(suite4)

unittest.TextTestRunner(verbosity=0).run(suite)


