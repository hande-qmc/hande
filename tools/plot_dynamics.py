#!/usr/bin/env python
'''plot_dynamics.py file

Plot the population dynamics and energy profile of a HANDE QMC output file.'''

import sys
import pandas as pd
import numpy as np
import pyhande
import matplotlib.pyplot as pyplot

def main(f):
    # [review] - JSS: either have or don't have a docstring---don't just add a null one!
    '''
''' 
    out = pyhande.extract.extract_data(f)
    data = out[1]
    shoulder = pyhande.analysis.shoulder_estimator(data)

    # Plot the total poulation over the entire range
    pyplot.subplot(3,1,1)
    pyplot.plot(data['iterations'], data['# H psips'])
    # [review] - JSS: spelling
    pyplot.xlabel('itteration')
    pyplot.ylabel('Total Population')

    # Plot the energy estimators over the entire range
    pyplot.subplot(3,1,2)
    pyplot.plot(data['iterations'], data['\sum H_0j N_j']/data['N_0'], label='Proj. Energy')
    pyplot.plot(data['iterations'], data['Shift'], label='Shift')
    # [review] - JSS: spelling
    pyplot.xlabel('itteration')
    pyplot.ylabel('Energy / $E_{h}$')
    pyplot.legend()

    pyplot.subplot(3,1,3)
    # Plot the total population up to the shoulder
    height = shoulder['mean']['shoulder height']
    data_around_shoulder = data[np.logical_and(data['# H psips'] < 1.1*height, 
                           data['# H psips'] > 0.9*height) ]
    pyplot.plot(data_around_shoulder['iterations'], 
           data_around_shoulder['# H psips'], label='Total Population')
    x_points = [min(data_around_shoulder['iterations']), 
                max(data_around_shoulder['iterations'])]
    pyplot.plot(x_points, [height, height], label='Shoulder Height') 
    # [review] - JSS: spelling
    pyplot.xlabel('itteration')
    pyplot.ylabel('Population')
    pyplot.legend(loc=2)
    pyplot.draw()
    pyplot.show()

    # Also print out the infomation about the shoulder
    # Stealing from reblock_hande.py
    try:
        float_fmt = '{0:-#.8e}'.format
        float_fmt(1.0)
    except ValueError:
        # GAH.  Alternate formatting only added to format function after
        # python 2.6..
        float_fmt = '{0:-.8e}'.format
    print(shoulder.to_string(float_format=float_fmt, line_width=80))

if __name__ == '__main__':

    # [review] - JSS: At the very least handle the case where a filename is not supplied or the
    # [review] - JSS: first argument is --help or -h.  I would use argparse for this (see
    # [review] - JSS: tests/extract_test_data.py).
    main(sys.argv[1])
