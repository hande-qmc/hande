#!/usr/bin/python
'''average_beta.py [options] file1 file2 ... fileN

Average over different DMQMC calculations contained in the output files, file1 file2 ... fileN, at each inverse temperature (beta) contained in the output files.'''

import optparse
import re
import scipy
import scipy.stats
import sys
try:
    import matplotlib.pyplot as pyplot
    USE_MATPLOTLIB = True
except ImportError:
    USE_MATPLOTLIB = False

# data columns
BETA_COL = 0
SHIFT_COL = 1
TR_HRHO_COL = 4
TR_RHO_COL = 3

class Stats:
    def __init__(self, mean, se):
        self.mean = mean
        self.se = se

class Data:
    def __init__(self, shift, Tr_Hrho, Tr_rho, init_lists=True):
        if init_lists:
            self.shift = [shift]
            self.Tr_Hrho = [Tr_Hrho]
            self.Tr_rho = [Tr_rho]
        else:
            self.shift = shift
            self.Tr_Hrho = Tr_Hrho
            self.Tr_rho = Tr_rho
    def append(self, shift, Tr_Hrho, Tr_rho):
        self.shift.append(shift)
        self.Tr_Hrho.append(Tr_Hrho)
        self.Tr_rho.append(Tr_rho)
    def calculate_stats(self):
        s = Data(0,0,0,False)
        for attr in ('shift', 'Tr_Hrho', 'Tr_rho'):
            mean = scipy.mean(getattr(self, attr))
            se = scipy.stats.sem(getattr(self, attr))
            setattr(s, attr, Stats(mean, se))
        return s
        
def stats_ratio(s1, s2):
    mean = s1.mean/s2.mean
    se = scipy.sqrt( (s1.se/s1.mean)**2 + (s2.se/s2.mean)**2 )*mean
    return Stats(mean, se)

def extract_data(data_files):

    data = {}
    
    # ignore lines beginning with # as the first non-whitespace character.
    comment = '^ *#'
    # start of data table
    start_regex = '^ # iterations'
    # end of data table is a blank line
    end_regex = '^ *$'
    # timestep value from echoed input data
    tau_regex = r'(?<=\btau\b)(.*)'

    for data_file in data_files:

        f = open(data_file)

        have_data = False
        for line in f:

            if re.search(end_regex, line):
                have_data = False
            elif have_data:
                # extract data!
                # data structure:
                #   data[beta] = ( [data from runs], [data item from runs], ...., )
                if not re.match(comment, line):
                    words = line.split()
                    beta = tau*float(words[BETA_COL])
                    shift = float(words[SHIFT_COL])
                    Tr_Hrho = float(words[TR_HRHO_COL])
                    Tr_rho = float(words[TR_RHO_COL])
                    if beta in data:
                        data[beta].append(shift, Tr_Hrho, Tr_rho)
                    else:
                        data[beta] = Data(shift, Tr_Hrho, Tr_rho)
            elif re.search(start_regex, line):
                have_data = True
            elif re.search(tau_regex, line):
                # get tau
                tau = float(re.search(tau_regex, line).group())

        f.close()

    return data

def get_data_stats(data):

    stats = {}
    for (beta,data_values) in data.iteritems():

        stats[beta] = data_values.calculate_stats()
        # Bonus!  Now calculate ratio Tr[H\rho]/Tr[\rho]
        stats[beta].E = stats_ratio(stats[beta].Tr_Hrho, stats[beta].Tr_rho)

    return stats

def print_stats(stats, trace=True, shift=False):

    print r'#     \beta     ',
    if shift:
        print '      shift           shift s.e.',
    if trace:
        print '      Tr[Hp]/Tr[p]    Tr[Hp]/Tr[p] s.e.',
    print
    for beta in sorted(stats.iterkeys()):

        data = stats[beta]
        print '%16.8f' % (beta) , 
        if shift:
            print '%16.8f%16.8f' % (data.shift.mean, data.shift.se) ,
        if trace:
            print '%16.8f%16.8f' % (data.E.mean, data.E.se) ,
        print

def plot_stats(stats, trace=True, shift=False):

    if USE_MATPLOTLIB:

        beta = sorted(stats.keys())
        if shift:
            # data
            data = [stats[b].shift.mean for b in beta] 
            # errorbars
            err = [stats[b].shift.se for b in beta] 
            pyplot.errorbar(beta, data, yerr=err, label='S')
        if trace:
            # data
            data = [stats[b].E.mean for b in beta] 
            # errorbars
            err = [stats[b].E.se for b in beta] 
            pyplot.errorbar(beta, data, yerr=err, label='Tr[Hp]/Tr[p]')
            pyplot.legend()
            pyplot.show()

def parse_options(args):

    parser = optparse.OptionParser(usage = __doc__)
    parser.add_option('--plot', action='store_true', default=False, help='Plot averaged data.')
    parser.add_option('--with-shift', action='store_true', dest='with_shift', default=False, help='Analyse shift data.')
    parser.add_option('--without-shift', action='store_false', dest='with_shift', help='Do not analyse shift data.  Default.')
    parser.add_option('--with-trace', action='store_true', dest='with_trace', default=True, help='Analyse trace data.  Default.')
    parser.add_option('--without-trace', action='store_false', dest='with_trace', help='Do not analyse trace data.')

    (options, filenames) = parser.parse_args(args)
    
    if len(filenames) == 0:
        parser.print_help()
        sys.exit(1)

    return (options, filenames)

if __name__ == '__main__':

    (options, data_files) = parse_options(sys.argv[1:])
    data = extract_data(data_files)
    stats = get_data_stats(data)
    print_stats(stats, options.with_trace, options.with_shift)
    if options.plot:
        plot_stats(stats, options.with_trace, options.with_shift)
