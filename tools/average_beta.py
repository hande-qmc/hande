#!/usr/bin/python
'''average_beta.py [options] file1 file2 ... fileN

Average over different DMQMC calculations contained in the output files, file1 file2 ... fileN, at each inverse temperature (beta) contained in the output files.'''

import optparse
import re
import scipy
import scipy.stats
import scipy.interpolate
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
    def __init__(self, shift, Tr_rho, numerators, init_lists=True):
        if init_lists:
            self.shift = [shift]
            self.Tr_rho = [Tr_rho]
            self.numerators = [numerators]
        else:
            self.shift = shift
            self.Tr_rho = Tr_rho
            self.numerators = numerators
    def append(self, shift, Tr_rho, numerators):
        self.shift.append(shift)
        self.Tr_rho.append(Tr_rho)
        self.numerators.append(numerators)
    def calculate_stats(self):
        s = Data(0,0,0,False)
        for attr in ('shift', 'Tr_rho'):
            mean = scipy.mean(getattr(self, attr))
            se = scipy.stats.sem(getattr(self, attr))
            setattr(s, attr, Stats(mean, se))
        mean = scipy.mean(self.numerators,0) 
        se = scipy.stats.sem(self.numerators,0)
        setattr(s, 'numerators', Stats(mean,se))

        return s

def get_estimator_headings(data_files):
    # start of data table
    start_regex = '^ # iterations'
    data_file = data_files[0]
    n_numerators = calc_num_numerators(data_file)
    in_file = open(data_file)
    estimator_headings = []
    while in_file:
       line = in_file.next()
       if re.search(start_regex, line): 
           headings = line.split()
           for i in range(7,7+n_numerators):
              estimator_headings.append(headings[i])
           for j in range(0,len(estimator_headings)):
              estimator_headings[j] = 'Tr['+estimator_headings[j][13]+'p]/Tr[p]'
           break
           
    return estimator_headings

def calc_num_numerators(data_file):
    # start of data table
    start_regex = '^ # iterations'
    in_file = open(data_file)
    n_numerators = 0
    while in_file:
       line = in_file.next()
       if re.search(start_regex, line):
           line = in_file.next()   
           n_columns = len(line.split())
           n_numerators = n_columns - 7
           break
    return n_numerators
        
def stats_ratio(s1, s2):
    mean = s1.mean/s2.mean
    se = scipy.sqrt( (s1.se/s1.mean)**2 + (s2.se/s2.mean)**2 )*mean
    return Stats(mean, se)

def stats_array_ratio(numerator_array,denominator):
    mean = []
    se = []
    for i in range(0,numerator_array.mean.size):
       mean.append(numerator_array.mean[i]/denominator.mean)
       se.append(scipy.sqrt( (numerator_array.se[i]/numerator_array.mean[i])**2 + (denominator.se/denominator.mean)**2 )*abs(mean[i]))
        
    return Stats(mean,se)

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
       
        num_numerators = calc_num_numerators(data_file)
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
                    numerators = []
                    for i in range(4,4+num_numerators):
                       numerators.append(float(words[i])) 
                    
                    Tr_rho = float(words[TR_RHO_COL])
                    if beta in data:
                        data[beta].append(shift, Tr_rho, numerators)
                    else:
                        data[beta] = Data(shift, Tr_rho, numerators)
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
        stats[beta].estimators = stats_array_ratio(stats[beta].numerators, stats[beta].Tr_rho)
    return stats

#def calculate_energy_fit(energies, betas, order):
#    poly_coeffs = scipy.polyfit(betas, energies, order)
#    energy_fit = scipy.polyval(poly_coeffs,betas)
#    return energy_fit

def calculate_spline_fit(energies, betas, weights):
    smoothing_condition = len(energies)-scipy.sqrt(2*len(energies))
    tck = scipy.interpolate.splrep(betas, energies, weights, k=3, s=smoothing_condition)
    
    spline_fit = scipy.interpolate.splev(betas,tck)
    spline_fit_deriv = scipy.interpolate.splev(betas,tck,der=1)
    
    return spline_fit, spline_fit_deriv

def calculate_specific_heat(energies, betas, weights):
#    poly_coeffs = scipy.polyfit(betas, energies, order)
#    diff_poly_coeffs = differentiate_polynomial(poly_coeffs)
#    specific_heat = scipy.polyval(diff_poly_coeffs,betas)
    energy_fit, specific_heat = calculate_spline_fit(energies, betas, weights)
    for i in range(0,len(specific_heat)):
        specific_heat[i] = -1.*betas[i]*betas[i]*specific_heat[i]
    return energy_fit, specific_heat

def print_stats(stats, estimator_headings, trace=True, shift=False, heat_capacity=False):
    energies = []
    betas = []
    weights = []
    if heat_capacity:
        for beta in sorted(stats.iterkeys()):
            data = stats[beta]
            weights.append(1./data.estimators.se[0])
            betas.append(beta)
            energies.append(data.estimators.mean[0]) 
        energy_fit, specific_heats = calculate_specific_heat(energies, betas, weights)
    print '#           beta    ',
    if shift:
        print 'shift           shift s.e.    ',
    if trace:
        for i in range(0,len(estimator_headings)):
            print estimator_headings[i]+'            s.e.    ',
    if heat_capacity:    
        print '  Energy Fit   Heat Capacity'
    print
    counter = 0
    for beta in sorted(stats.iterkeys()):

        data = stats[beta]
        print '%16.8f' % (beta) ,
        if shift:
            print '%16.8f%16.8f' % (data.shift.mean, data.shift.se) ,
        if trace:
            for i in range(0,len(data.estimators.mean)):
                print '%16.8f%16.8f' % (data.estimators.mean[i], data.estimators.se[i]) ,
        if heat_capacity:
            print '%16.8f%16.8f' % (energy_fit[counter], specific_heats[counter]) ,
        counter = counter + 1
        print

#def differentiate_polynomial(poly_coeffs):
#    diff_poly_coeffs = []
#    for i in range(0,len(poly_coeffs)-1):   
#        diff_poly_coeffs.append(poly_coeffs[i]*(len(poly_coeffs)-i-1))
#    return diff_poly_coeffs


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
            data = [stats[b].estimators.mean[0] for b in beta] 
            # errorbars
            err = [stats[b].estimators.se[0] for b in beta] 
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
    parser.add_option('--with-heat_capacity', action='store_true', dest='with_heat_capacity', default=False, help='Calculate heat capacity')
    parser.add_option('--without-heat_capacity', action='store_false', dest='with_heat_capacity', help='Do not calcualate heat capacity. Default')
    (options, filenames) = parser.parse_args(args)
    
    if len(filenames) == 0:
        parser.print_help()
        sys.exit(1)

    return (options, filenames)

if __name__ == '__main__':

    (options, data_files) = parse_options(sys.argv[1:])
    estimator_headings = get_estimator_headings(data_files)
    print estimator_headings
    data = extract_data(data_files)
    stats = get_data_stats(data)
    print_stats(stats, estimator_headings, options.with_trace, options.with_shift, options.with_heat_capacity)
    if options.plot:
        plot_stats(stats, options.with_trace, options.with_shift)
