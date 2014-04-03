import unittest

import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
import pyblock

class BlockingTests(unittest.TestCase):
    pass

def main():
    unittest.main()

if __name__ == '__main__':

    main()
