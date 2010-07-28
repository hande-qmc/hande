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
se: the standard error of the data (=std. dev./sqrt(block_size));
se_error: the estimate of the error associated with the standard error.
'''
    def __init__(self, block_size, mean, se, se_error):

        self.block_size = block_size
        self.mean = mean
        self.se = se
        self.se_error = se_error
    def __repr__(self):

        return (self.block_size, self.mean, self.se, self.se_error).__repr__()


class Data(object):
    '''Class for holding the data during the blocking process.

data_col: index (starting from 0) of the column containing the data.'''

    def __init__(self, data_col):

        self.data_col = data_col

        # Initialise some null values we'll use during the analysis.
        self.data = []
        self.stats = []

    def add_to_data(self, value):
        '''Add value to the list of data.'''

        self.data.append(value)

    def sort_by_index(self, indices):
        '''Sort data according to the (unsorted) list of indices.'''

        combined = zip(indices, self.data)
        combined.sort()
        (sorted_indices, self.data) = zip(*combined)

    def reblock(self):
        '''Reblock the data by successively taking the mean of adjacent data points.'''

        block_size = len(self.data)/2
        self.data = [0.5*(self.data[2*i]+self.data[2*i+1]) for i in range(block_size)]

    def add_stats(self):
        '''Calculate the statistics of the current set of data and append to the stats list.'''

        self.stats.append(self.calculate_stats())

    def calculate_stats(self):
        '''Calculate the statistics of the current set of data and return the corresponding Stats object.'''

        mean = self.calculate_mean()
        (se, se_err) = self.calculate_se(mean)
        block_size = len(self.data)
        return Stats(block_size, mean, se, se_err)

    def calculate_mean(self):
        '''Calculate the mean of the current set of data.'''

        sum = 0
        for value in self.data:
            sum += value
        return sum/len(self.data)

    def calculate_se(self, mean):
        '''Calculate the standard error and estimate associated error of the current set of data.'''

        # See "Error estimates on averages of correlated data" by Flyvbjerg and Petersen, JCP 91 461 (1989).
        # They give the standard error and associated error in the standard error to be:
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
        se = sqrt(c0/(size*(size-1)))
        se_err = se*1.0/(sqrt(2*(size-1)))
        return (se, se_err)


class DataBlocker(object):
    '''Class for performing a blocking analysis.

datafiles: list of datafiles containing the data to be extracted and blocked;
start_regex: regular expression indicating that subsequent lines contain the data;
end_regex: regular expression indicating that the end of a block of data
           (multiple data blocks can occur within one file);
index_col: index (starting from 0) of the column containing the index of the
           data (e.g. Monte Carlo cycle) to be blocked.  The index is only used to
           initially sort the data;
data_cols: list of indices (starting from 0) of the columns containing the data;
start_index: block only data with an index greater than or equal to start_index.
block_all: assume all lines contain data apart from comment lines.  The regular expressions are ignored if this is true.
'''
    def __init__(self, datafiles, start_regex, end_regex, index_col, data_cols, start_index, block_all=False, combination=False):

        self.datafiles = datafiles
        self.start_regex = re.compile(start_regex)
        self.end_regex = re.compile(end_regex)
        self.comment_regex = re.compile('^ *#')
        self.index_col = index_col
        self.start_index = start_index
        self.block_all = block_all

        self.data = [Data(data_col) for data_col in data_cols]

        self.combination = combination
        if self.combination == '/':
            self.combination_fn = self.calculate_combination_division
        elif self.combination:
            raise Exception, 'Unimplemented combination: %s.' % (self.combination)

        # Dictionaries for storing the covariance and combination stats between
        # two different data sets at each given block size.
        # The key i,j is used to distinguish different combinations, where
        # i and j correspond to values in data_cols.
        # These are symmetric and thus we choose to work using i<j.
        self.covariance = {}
        self.combination_stats = {}

    def get_data(self):
        '''Extract the relevant data from the datafiles.'''

        indices = []
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
                    index = float(d[self.index_col])
                    if index >= self.start_index:
                        indices.append(index)
                        for data in self.data:
                            data.add_to_data(float(d[data.data_col]))
                # are we about to start receiving the data?
                if re.match(self.start_regex, line):
                    have_data = True
        # Now we sort the data according to the index.
        for data in self.data:
            data.sort_by_index(indices)

    def blocking(self):
        '''Perform a blocking analysis on the data.
        
This destroys the data stored in self.data.data'''

        # See "Error estimates on averages of correlated data" by Flyvbjerg and Petersen, JCP 91 461 (1989).
        block_size = len(self.data[0].data)
        while block_size >= 2:

            for (i, data) in enumerate(self.data):
                # calculate stats for this data set
                data.add_stats()

            if len(self.data) > 1:
                # Bonus: also calculate the covariance.
                for i in range(len(self.data)):
                    for j in range(i+1, len(self.data)):
                        key = '%s,%s' % (self.data[i].data_col, self.data[j].data_col)
                        if key not in self.covariance:
                            self.covariance[key] = []
                        self.covariance[key].append(self.calculate_covariance(i, j))
                        # Added bonus: calculate combination if desired.
                        if self.combination == '/':
                            if key not in self.combination_stats:
                                self.combination_stats[key] = []
                            self.combination_stats[key].append(self.calculate_combination_division(i,j))

            for (i, data) in enumerate(self.data):
                # Update length of block size after this reblocking cycle.
                if i == 0:
                    block_size /= 2

                data.reblock()

    def calculate_covariance(self, i, j):
        '''Calculate the covariance between the i-th data item and the j-th data item.
        
Note that this assumes that the mean stored in the Data class corresponds to
that of the current set of data.'''

        # cov(X,Y) = E[(X-\mu_X)(Y-\mu_Y)]
        #          = 1/N \sum_i=1^N (X_i - \mu_X)(Y_i -\mu_Y)

        cov = 0
        for x in range(len(self.data[i].data)):
            cov += (self.data[i].data[x] - self.data[i].stats[-1].mean)*(self.data[j].data[x] - self.data[j].stats[-1].mean)
        return cov/len(self.data[i].data)

    def calculate_combination_division(self, i, j):
        '''Find the mean and standard error of f, where f = X_i/X_j, where X_i is the i-th data set and similarly for X_j.

i,j: elements in the data array corresponding to X_i and X_j.

NOTE: we assume that the statistical values in data[i].stats are consistent
with the data in data[i].data and similarly for data[j].data and the covaraince
data (ie this must be called after add_stats and calculate_covariance for each
blocking cycle).

The mean of f, <f>, can simply be found from <X_i>/<X_j>.

Denoting the standard deviation of the data set X as s(X) and the standard error of data set X as se(X),
the standard deviation of f can be evaluated using:

    (s(f)/f)^2 = (s(X_i)/<X_i>)^2 + (s(X_i)/<X_i>)^2 - 2 cov(X_i,X_j)/(<X_i> <X_j>)

where cov(X_i,X_j) is the covariance between the two data sets.  The standard error then follows:
    
    se(f) = s(f)/\sqrt(N)
    
where N is the number of elements in the data set and hence we obtain:
    
    se(f) = f [ (se(X_i)/<X_i>)^2 + (se(X_j)/<X_j>)^2 - 2 cov(X_i,X_j)/(N <X_i> <X_j>) ]^1/2.'''

        meani = self.data[i].stats[-1].mean
        meanj = self.data[j].stats[-1].mean
        sei = self.data[i].stats[-1].se
        sej = self.data[j].stats[-1].se
        cov_key = '%s,%s' % tuple(sorted((self.data[i].data_col,self.data[j].data_col)))
        cov = self.covariance[cov_key][-1]
        nblocks = self.data[j].stats[-1].block_size

        meanf = meani/meanj

        sef = abs(meanf*sqrt( (sei/meani)**2 + (sej/meanj)**2 - 2*cov/(nblocks*meani*meanj) ))

        return Stats(nblocks, meanf, sef, 0)

    def show_blocking(self, plotfile=''):
        '''Print out the blocking data and show a graph of the behaviour of the standard error with block size.
        
If plotfile is given, then the graph is saved to the specifed file rather than being shown on screen.'''

        # print blocking output
        # header...
        print '%-11s' % ('# of blocks'),
        fmt = '%-14s %-12s %-18s '
        header = ('mean (X_%i)', 'std.err. (X_%i)', 'std.err.err. (X_%i)')
        for data in self.data:
            data_header = tuple(x % (data.data_col) for x in header)
            print fmt % data_header,
        for key in self.covariance:
            str = 'cov(X_%s,X_%s)' % tuple(key.split(','))
            print '%-12s' % (str),
        for key in self.combination_stats:
            fmt = ['mean (X_%s'+self.combination+'X_%s)', 'std.err. (X_%s'+self.combination+'X_%s)']
            strs = tuple([s % tuple(key.split(',')) for s in fmt])
            print '%-16s %-18s' % strs,
        print
        # data
        block_fmt = '%-11i'
        fmt = '%-#14.12g %-#12.8e %-#18.8e '
        for s in range(len(self.data[0].stats)):
            print block_fmt % (self.data[0].stats[s].block_size),
            for data in self.data:
                print fmt % (data.stats[s].mean, data.stats[s].se, data.stats[s].se_error),
            for cov in self.covariance.itervalues():
                print '%-#12.9g' % (cov[s]),
            for comb in self.combination_stats.itervalues():
                print '%-#16.12g %-#18.12g' % (comb[s].mean, comb[s].se),
            print

        # plot standard error 
        if PYLAB:
            # one sub plot per data set.
            nplots = len(self.data)
            for (i, data) in enumerate(self.data):
                pylab.subplot(nplots, 1, i+1)
                blocks = [stat.block_size for stat in data.stats]
                se = [stat.se for stat in data.stats]
                se_error = [stat.se_error for stat in data.stats]
                pylab.semilogx(blocks, se, 'g-', basex=2, label=r'$\sigma(X_%s)$' % (data.data_col))
                pylab.errorbar(blocks, se, yerr=se_error, fmt=None, ecolor='g')
                xmax = 2**pylab.ceil(pylab.log2(blocks[0]+1))
                pylab.xlim(xmax, 1)
                pylab.ylabel('Standard error')
                pylab.legend(loc=2)
                if i != nplots - 1:
                    # Don't label x axis points.
                    ax = pylab.gca()
                    ax.set_xticklabels([])
            pylab.xlabel('# of blocks')
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
    parser.add_option('-d', '--data', dest='data_cols', type='int', default=[], action='append', help='Set the column(s) (starting from 0) containing the data items. If more than one column is given, the covariance between the sets of data is also calculated.  Default: 1.')
    parser.add_option('-f', '--from', dest='start_index', type='int', default=0, help='Set the index from which the data is blocked.  Data with a smaller index is discarded.  Default: %default.')
    parser.add_option('-p', '--plotfile', help='Save a plot of the blocking analysis to PLOTFILE rather than showing the plot on screen (default behaviour).')
    parser.add_option('-o','--operation', help='')

    (options, filenames) = parser.parse_args(args)

    # Set additional defaults.
    if not options.data_cols:
        options.data_cols = [1]

    if len(filenames) == 0:
        parser.print_help()
        sys.exit(1)

    return (options, filenames)

if __name__ == '__main__':
    (options, filenames) = parse_options(sys.argv[1:])

    my_data = DataBlocker(filenames, options.start_regex, options.end_regex, options.index_col, options.data_cols, options.start_index, options.all, options.operation)

    my_data.get_data()
    my_data.blocking()
    my_data.show_blocking(options.plotfile)
