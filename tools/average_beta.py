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
TR_H2RHO_COL = 5
TR_RHO_COL = 3
H2_is_present = False
H_is_present = False

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
    def calculate_covs(self):
        covs = []
        for i in range(0,len(self.numerators[0])):
            covariance_matrix = scipy.cov(scipy.vstack((scipy.array(self.numerators)[:,i],self.Tr_rho)))
            covs.append(covariance_matrix[0][1]/len(self.numerators))
        return covs

def get_estimator_headings(data_files):
# Get collumn headings from input file
    start_regex = '^ # iterations'
    data_file = data_files[0]
    in_file = open(data_file)
    estimator_headings = []
    estimator_col_no = []
    while in_file:
       line = in_file.next()
       if re.search(start_regex, line): 
           headings = line.split()
           index = 0
           for k in range(0,len(headings)):
               if headings[k] == '\\sum\\rho_{ij}H_{ji}':
                   estimator_col_no.append(k-3)
                   estimator_headings.append('Tr[Hp]/Tr[p]')
                   H_is_present = True
                   TR_HRHO_INDEX = index
                   index = index + 1
               elif headings[k] == '\\sum\\rho_{ij}S_{ji}':
                   estimator_col_no.append(k-3)
                   estimator_headings.append('Tr[Sp]/Tr[p]')
                   index = index + 1
               elif headings[k] == '\\sum\\rho_{ij}C_{ji}':
                   estimator_col_no.append(k-3)
                   estimator_headings.append('Tr[Cp]/Tr[p]')
                   index = index + 1 
               elif headings[k] == '\\sum\\rho_{ij}H2_{ji}':
                   estimator_col_no.append(k-3)
                   estimator_headings.append('Tr[H2p]/Tr[p]')
                   H2_is_present = True
                   TR_H2RHO_INDEX = index
                   index = index + 1
               else:
                   print '# Warning: Column name, '+headings[k]+', not recognised...'
           break
           
    return estimator_headings, estimator_col_no

def stats_ratio(s1, s2):
    mean = s1.mean/s2.mean
    se = scipy.sqrt( (s1.se/s1.mean)**2 + (s2.se/s2.mean)**2 )*abs(mean)
    return Stats(mean, se)

def stats_array_ratio(numerator_array,denominator, covs):
# Calculate the new stats for numerators divided by tr(rho)
    mean = []
    se = []
    # Loop over numerators and calculate s.e and mean for this estimator
    for i in range(0,numerator_array.mean.size):
       mean.append(numerator_array.mean[i]/denominator.mean)
       se.append(scipy.sqrt( (numerator_array.se[i]/numerator_array.mean[i])**2 + (denominator.se/denominator.mean)**2 - (2*covs[i]/(numerator_array.mean[i]*denominator.mean)) )*abs(mean[i]))
    return Stats(mean,se)

def stats_stochastic_specific_heat(beta, estimator_array):
    stochastic_specific_heats = []
    mean = (beta**2)*(estimator_array.mean[TR_H2RHO_INDEX]-estimator_array.mean[TR_HRHO_INDEX]**2)
    se = scipy.sqrt((beta**4)*(estimator_array.se[TR_H2RHO_INDEX]**2+4*(estimator_array.mean[TR_HRHO_INDEX]**2)*estimator_array.se[TR_HRHO_INDEX]**2))
    return Stats(mean,se)

def extract_data(data_files, estimtor_col_no):
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
                    numerators = []
                    for i in range(0,len(estimtor_col_no)):
                       numerators.append(float(words[estimtor_col_no[i]])) 
                    
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
        covariances = data_values.calculate_covs()
        stats[beta] = data_values.calculate_stats()
        
        # Bonus!  Now calculate ratio Tr[H\rho]/Tr[\rho]
        stats[beta].estimators = stats_array_ratio(stats[beta].numerators, stats[beta].Tr_rho, covariances)
        if H2_is_present and H_is_present:
            stats[beta].stochastic_specific_heat = stats_stochastic_specific_heat(beta,stats[beta].estimators)
  
    return stats

def calculate_spline_fit(energies, betas, weights):
# Perform a cubic spline best fit on the energy data. 
# weights: array of 1/(the standard deviations for each E(beta)), Used for calculating
#          the least-squares spline fit
    
    # Choose default smoothing parameter m-sqrt(2*m) < s < m+sqrt(2*m). 
    # Where m is number of data points. 
    # Defaults is s = m - sqrt(2*m)

    # Caluclate knots, spline co-efficients and degree of spline (t,c,k=3)
    tck = scipy.interpolate.splrep(betas, energies, weights, k=3)
    # Calculate the spline fit using (t,c,k)
    spline_fit = scipy.interpolate.splev(betas,tck)
    # Calculate the first derivateive of spline fit using (t,c,k) 
    spline_fit_deriv = scipy.interpolate.splev(betas,tck,der=1)
    
    return spline_fit, spline_fit_deriv

def calculate_specific_heat(energies, betas, weights):
    # Calculate energy spline fit and first derivative
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
    if H2_is_present and H_is_present:
       print 'Stoch. Spec. Heat  s.e.    '
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
        if H2_is_present and H_is_present:
            print '%16.8f%16.8f' % (data.stochastic_specific_heat.mean, stochastic_specific_heat.se) ,
        if heat_capacity:
            print '%16.8f%16.8f' % (energy_fit[counter], specific_heats[counter]) ,
        counter = counter + 1
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
    estimator_headings, estimtor_col_no = get_estimator_headings(data_files)
    data = extract_data(data_files, estimtor_col_no)
    stats = get_data_stats(data)
    print_stats(stats, estimator_headings, options.with_trace, options.with_shift, options.with_heat_capacity)
    if options.plot:
        plot_stats(stats, options.with_trace, options.with_shift)
