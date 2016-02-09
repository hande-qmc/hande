#!/usr/bin/env python

import os
import sys
try:
    import pyhande
except ImportError:
    # Try to find pyhande in the hande directory tree.
    _script_dir = os.path.dirname(os.path.abspath(__file__))
    sys.path.append(os.path.join(_script_dir, '../pyhande'))
    import pyhande
import pyhande.testcode


def main(filename):

    data = pyhande.testcode.extract_test_data(filename)
    output = []
    for (key, val) in data.items():
        output.append((key, val.to_string(index=False, index_names=False, float_format='{:.10f}'.format)))

    return '\n\n'.join('\n'.join(calc_out) for calc_out in output)

if __name__ == '__main__':

    if len(sys.argv) == 1:
        print('usage: extract_test_data.py file')
        sys.exit(1)
    print(main(sys.argv[1]))
