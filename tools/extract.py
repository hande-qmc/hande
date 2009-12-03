#!/usr/bin/env python

import re
import sys

def extract(filename):

    f = open(filename)
    data = f.readlines()
    f.close()

    for line in data:
        if re.search('Ground state:', line):
            print "ground_state\n%s" % (line.split()[-1])

if __name__ == '__main__':
    filename = sys.argv[1]
    extract(filename)
