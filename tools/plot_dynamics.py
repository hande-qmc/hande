#!/usr/bin/env python
'''plot_dynamics.py file

Plot the population dynamics and energy profile of a HANDE QMC output file.'''

import sys
import pandas as pd
import numpy as np
import pyhande
import matplotlib.pyplot as pyplot
import argparse

def main(f):

    out = pyhande.extract.extract_data(f)
    data = out[1]
    shoulder = pyhande.analysis.plateau_estimator(data)

    # Plot the total population over the entire range
    pyplot.subplot(3,1,1)
    pyplot.plot(data['iterations'], data['# H psips'])
    # [review] - JSS: spelling
    pyplot.xlabel('iteration')
    pyplot.ylabel('Total Population')

    # Plot the energy estimators over the entire range
    pyplot.subplot(3,1,2)
    pyplot.plot(data['iterations'], data['\sum H_0j N_j']/data['N_0'], label='Proj. Energy')
    pyplot.plot(data['iterations'], data['Shift'], label='Shift')
    # [review] - JSS: spelling
    pyplot.xlabel('iteration')
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
    pyplot.xlabel('iteration')
    pyplot.ylabel('Population')
    pyplot.legend(loc=2)
    pyplot.draw()
    pyplot.show()

    # Also print out the information about the shoulder
    # Stealing from reblock_hande.py
    try:
        float_fmt = '{0:-#.8e}'.format
        float_fmt(1.0)
    except ValueError:
        # GAH.  Alternate formatting only added to format function after
        # python 2.6..
        float_fmt = '{0:-.8e}'.format
    print(shoulder.to_string(float_format=float_fmt, line_width=80))

def parse_args(args):

    parser = argparse.ArgumentParser(description='Plots the population and energy estimators of an FCIQMC/CCMC calulation')
    parser.add_argument('-f', '--file', action='store', type=str, help='File to plot.')
    opts = parser.parse_args(args)
    if not opts.file:
        parser.print_usage()
        print('Must supply a file to plot')
        sys.exit()
    return opts.file

if __name__ == '__main__':

    file = parse_args(sys.argv[1:])
    main(file)
