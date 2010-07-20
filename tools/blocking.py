#!/usr/bin/python
'''

blocking.py [options] file1 file2 ... fileN

Perform blocking analysis, where file1 file2 ... fileN are files containing the
data to be analysed.  An arbitrary (positive) number of files can be analysed.
Lines with a # as the first non-space character are treated as comments and are ignored.'''

from math import sqrt
import operator
import optparse
import re
import sys
try:
    import pylab
    PYLAB = True
except ImportError:
    print "Can't import matplotlib.  Skipping graph production."
    PYLAB = False

__author__ = 'James Spencer'

class Stats(object):
    '''Class to store statistics about a set of data.

block_size: the number of points in a set of data; 
mean: the mean of the data;
sd: the standard deviation of the data;
sd_error: the estimate of the error associated with the standard deviation.
'''
    def __init__(self, block_size, mean, sd, sd_error):

        self.block_size = block_size
        self.mean = mean
        self.sd = sd
        self.sd_error = sd_error
    def __repr__(self):

        return (self.block_size, self.mean, self.sd, self.sd_error).__repr__()

class DataBlocker(object):
    '''Class for performing a blocking analysis.
datafiles: list of datafiles containing the data to be extracted and blocked;
start_regex: regular expression indicating that subsequent lines contain the data;
end_regex: regular expression indicating that the end of a block of data
           (multiple data blocks can occur within one file);
index_col: index (starting from 0) of the column containing the index of the
           data (e.g. Monte Carlo cycle) to be blocked.  The index is only used to
           initially sort the data;
data_col: index (starting from 0) of the column containing the data;
start_index: block only data with an index greater than or equal to start_index.
block_all: assume all lines contain data apart from comment lines.  The regular expressions are ignored if this is true.
'''
    def __init__(self, datafiles, start_regex, end_regex, index_col, data_col, start_index, block_all=False):

        self.datafiles = datafiles
        self.start_regex = re.compile(start_regex)
        self.end_regex = re.compile(end_regex)
        self.comment_regex = re.compile('^ *#')
        self.index_col = index_col
        self.data_col = data_col
        self.start_index = start_index
        self.block_all = block_all

        # Initialise some null values we'll use during the analysis.
        self.data = []
        self.stats = []
        self.nblocks = 0

    def get_data(self):
        '''Extract the relevant data from the datafiles.'''

        data_items = []
        for file in self.datafiles:
            f = open(file, 'r')
            have_data = False
            for line in f:
                # have we hit the end of the data?
                if re.match(self.end_regex, line):
                    have_data = False
                # do we have data to extract?
                if (have_data or self.block_all) and not re.match(self.comment_regex, line):
                    d = line.split()
                    (index, value) = (float(d[self.index_col]), float(d[self.data_col]))
                    if index >= self.start_index:
                        data_items.append((index, value))
                # are we about to start receiving the data?
                if re.match(self.start_regex, line):
                    have_data = True
        # Now we sort the data according to the index.
        data_items.sort()
        # Don't care about the index of the item for blocking, just need the raw data.
        self.data = [d[1] for d in data_items]

    def calculate_stats(self):
        '''Calculate the statistics of the current set of data and return the corresponding Stats object.'''

        mean = self.calculate_mean()
        (sd, sd_err) = self.calculate_sd(mean)
        block_size = len(self.data)
        return Stats(block_size, mean, sd, sd_err)

    def calculate_mean(self):
        '''Calculate the mean of the current set of data.'''

        sum = 0
        for value in self.data:
            sum += value
        return sum/len(self.data)

    def calculate_sd(self, mean):
        '''Calculate the standard deviation and associated error of the current set of data.'''

        # See "Error estimates on averages of correlated data" by Flyvbjerg and Petersen, JCP 91 461 (1989).
        # They give the standard deviation and associated error in the standard deviation to be:
        #   \sigma \approx \sqrt{c0/n-1} ( 1 \pm 1/(\sqrt(2(n-1))) )
        # where
        #   n is the number of items in the set of data;
        #   c0 is given by:
        #     c0 = 1/n \sum_k (x_k - \bar{x})^2
        # where
        #   x_k is the k-th data element;
        #   \bar{x} is the mean of the set of data.
        c0 = 0
        for value in self.data:
            c0 += (value - mean)**2
        size = len(self.data)
        sd = sqrt(c0/(size*(size-1)))
        sd_err = sd*1.0/(sqrt(2*(size-1)))
        return (sd, sd_err)

    def blocking(self):
        '''Perform a blocking analysis on the data.
        
This destroys the data stored in self.data'''

        # See "Error estimates on averages of correlated data" by Flyvbjerg and Petersen, JCP 91 461 (1989).
        block_size = len(self.data)
        while block_size >= 2:

            # calculate stats for this data set
            self.stats.append(self.calculate_stats())

            # reblock
            block_size /= 2
            self.data = [0.5*(self.data[2*i]+self.data[2*i+1]) for i in range(block_size)]

    def show_blocking(self, plotfile=''):
        '''Print out the blocking data and show a graph of the behaviour of the standard deviation with block size.  If plotfile is given, then the graph is saved to the specifed file rather than being shown on screen.'''

        # print blocking output
        fmt = '%-10s   %-16s  %-18s   %-24s'
        print fmt % ('block size', 'mean', 'standard deviation', 'standard deviation error')
        fmt = '%-10i  %-16.12g   %-18.12g   %-24.12g'
        for stat in self.stats:
            print fmt % (stat.block_size, stat.mean, stat.sd, stat.sd_error)

        # plot standard deviation
        if PYLAB:
            blocks = [stat.block_size for stat in self.stats]
            sd = [stat.sd for stat in self.stats]
            sd_error = [stat.sd_error for stat in self.stats]
            pylab.semilogx(blocks, sd, basex=2)
            pylab.errorbar(blocks, sd, yerr=sd_error)
            xmax = 2**pylab.ceil(pylab.log2(blocks[0]+1))
            pylab.xlim(xmax, 1)
            pylab.xlabel('Block size')
            pylab.ylabel('Standard deviation')
            if plotfile:
                pylab.savefig(plotfile)
            else:
                pylab.draw()
                pylab.show()

def parse_options(args):
    '''Parse command line options.'''

    parser = optparse.OptionParser(usage = __doc__)
    parser.add_option('-s', '--start', dest='start_regex', default='^ # iterations', help='Set the regular expression indicating the start of a data block.  Default: %default.')
    parser.add_option('-e', '--end', dest='end_regex', type='string', default=r'^ *$', help='Set the regular expression indicating the end of a data block.  Default: %default.')
    parser.add_option('-a', '--all', action='store_true', default=False, help='Assume all lines in the files contains data apart from comment lines. Regular expression options are ignored if --all is used.  Default: %default.')
    parser.add_option('-i', '--index', dest='index_col', type='int', default=0, help='Set the column (starting from 0) containing the index labelling each data item (e.g. number of Monte Carlo cycles). Default: %default.')
    parser.add_option('-d', '--data', dest='data_col', type='int', default=1, help='Set the column (starting from 0) containing the data items. Default: %default.')
    parser.add_option('-f', '--from', dest='start_index', type='int', default=0, help='Set the index from which the data is blocked.  Data with a smaller index is discarded.  Default: %default.')
    parser.add_option('-p', '--plotfile', help='Save a plot of the blocking analysis to PLOTFILE rather than showing the plot on screen (default behaviour).')

    (options, filenames) = parser.parse_args(args)

    if len(filenames) == 0:
        parser.print_help()
        sys.exit(1)

    return (options, filenames)

if __name__ == '__main__':
    (options, filenames) = parse_options(sys.argv[1:])

    my_data = DataBlocker(filenames, options.start_regex, options.end_regex, options.index_col, options.data_col, options.start_index, options.all)

    my_data.get_data()
    my_data.blocking()
    my_data.show_blocking(options.plotfile)
