#!/usr/bin/python
'''average_beta.py [options] file1 file2 ... fileN

Average over different DMQMC calculations contained in the output files, file1 file2 ... fileN, at each inverse temperature (beta) contained in the output files.'''

import optparse
import re
import scipy
import scipy.stats
import scipy.interpolate
import sys
from math import sqrt
from math import ceil
from math import floor
from math import exp
from math import log

try:
    import matplotlib.pyplot as pyplot
    USE_MATPLOTLIB = True
except ImportError:
    USE_MATPLOTLIB = False

# data columns
BETA_COL = 0
SHIFT_COL = 1
TR_RHO_COL = 2

TR_RHO2_is_present = False
H2_is_present = False
H_is_present = False
S_is_present = False
C_is_present = False
M2_is_present = False
Full_R2_is_present = False
RDM_R2_is_present = False

TR_RHO2_INDEX = -1
TR_SRHO_INDEX = -1
TR_CRHO_INDEX = -1
TR_M2RHO_INDEX = -1
TR_HRHO_INDEX = -1
TR_H2RHO_INDEX = -1
Full_R2_INDEX = -1
RDM_R2_INDEX = -1
RDM_TR1_INDEX = -1
RDM_TR2_INDEX = -1
PARTICLES_COL = 6


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
    global TR_RHO2_is_present
    global H_is_present
    global H2_is_present
    global S_is_present
    global C_is_present
    global M2_is_present
    global Full_R2_is_present
    global RDM_R2_is_present
    global TR_RHO2_INDEX
    global TR_HRHO_INDEX
    global TR_H2RHO_INDEX
    global TR_SRHO_INDEX
    global TR_CRHO_INDEX
    global TR_M2RHO_INDEX
    global Full_R2_INDEX
    global RDM_R2_INDEX
    global RDM_TR1_INDEX
    global RDM_TR2_INDEX
    global PARTICLES_COL
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
               if headings[k] == 'Trace_2':
                   estimator_col_no.append(k-TR_RHO_COL)
                   estimator_headings.append('Trace_2')
                   TR_RHO2_is_present = True
                   TR_RHO2_INDEX = index
                   index = index + 1
               elif headings[k] == '\\sum\\rho_{ij}H_{ji}':
                   estimator_col_no.append(k-TR_RHO_COL)
                   estimator_headings.append('Tr[Hp]/Tr[p]')
                   H_is_present = True
                   TR_HRHO_INDEX = index
                   index = index + 1
               elif headings[k] == '\\sum\\rho_{ij}S_{ji}':
                   estimator_col_no.append(k-TR_RHO_COL)
                   estimator_headings.append('Tr[Sp]/Tr[p]')
                   S_is_present = True
                   TR_SRHO_INDEX = index
                   index = index + 1
               elif headings[k] == '\\sum\\rho_{ij}C_{ji}':
                   estimator_col_no.append(k-TR_RHO_COL)
                   estimator_headings.append('Tr[Cp]/Tr[p]')
                   C_is_present = True
                   TR_CRHO_INDEX = index
                   index = index + 1 
               elif headings[k] == '\\sum\\rho_{ij}M2{ji}':
                   estimator_col_no.append(k-TR_RHO_COL)
                   estimator_headings.append('Tr[Mp]/Tr[p]')
                   M2_is_present = True
                   TR_M2RHO_INDEX = index
                   index = index + 1
               elif headings[k] == '\\sum\\rho_{ij}H2{ji}':
                   estimator_col_no.append(k-TR_RHO_COL)
                   estimator_headings.append('Tr[H2p]/Tr[p]')
                   H2_is_present = True
                   TR_H2RHO_INDEX = index
                   index = index + 1
               elif headings[k] == 'Full_R2_numerator':
                   estimator_col_no.append(k-TR_RHO_COL)
                   estimator_headings.append('Full_Renyi_2')
                   Full_R2_is_present = True
                   Full_R2_INDEX = index
                   index = index + 1
               elif headings[k] == 'Renyi_2_numerator_1':
                   estimator_col_no.append(k-TR_RHO_COL)
                   estimator_headings.append('RDM_Renyi_2')
                   RDM_R2_is_present = True
                   RDM_R2_INDEX = index
                   index = index + 1
               elif headings[k] == 'RDM1_trace_1':
                   estimator_col_no.append(k-TR_RHO_COL)
                   estimator_headings.append('RDM1_trace_1')
                   RDM_TR1_INDEX = index
                   index = index + 1
               elif headings[k] == 'RDM1_trace_2':
                   estimator_col_no.append(k-TR_RHO_COL)
                   estimator_headings.append('RDM1_trace_2')
                   RDM_TR2_INDEX = index
                   index = index + 1
           PARTICLES_COL = TR_RHO_COL + index + 1
           break
           
    return estimator_headings, estimator_col_no

def stats_ratio(s1, s2):
    mean = s1.mean/s2.mean
    se = scipy.sqrt( (s1.se/s1.mean)**2 + (s2.se/s2.mean)**2 )*abs(mean)
    return Stats(mean, se)

def stats_array_ratio(numerator_array, denominator, covs):
# Calculate the new stats for numerators divided by tr(rho)
    mean = []
    se = []
    # Loop over numerators and calculate s.e and mean for this estimator
    for i in range(0,numerator_array.mean.size):
       if i == RDM_R2_INDEX or i == RDM_TR1_INDEX or i == RDM_TR2_INDEX or i == TR_RHO2_INDEX or i == Full_R2_INDEX:
          mean.append(numerator_array.mean[i])
          se.append(numerator_array.se[i])
       else:
          mean.append(numerator_array.mean[i]/denominator.mean)
          se.append(scipy.sqrt( (numerator_array.se[i]/numerator_array.mean[i])**2 + (denominator.se/denominator.mean)**2 - (2*covs[i]/(numerator_array.mean[i]*denominator.mean)) )*abs(mean[i]))
    return Stats(mean,se)

def stats_stochastic_specific_heat(beta, estimator_array):
    mean = (beta**2)*(estimator_array.mean[TR_H2RHO_INDEX]-(estimator_array.mean[TR_HRHO_INDEX])**2)
    se = scipy.sqrt((beta**4)*(estimator_array.se[TR_H2RHO_INDEX]**2+4*(estimator_array.mean[TR_HRHO_INDEX]**2)*estimator_array.se[TR_HRHO_INDEX]**2))
    return Stats(mean,se)

def stats_full_renyi_2(estimator_array, trace_1):
    mean = -log(estimator_array.mean[Full_R2_INDEX]/(trace_1.mean*estimator_array.mean[TR_RHO2_INDEX]))/log(2)
    se = (1/log(2))*scipy.sqrt((estimator_array.se[Full_R2_INDEX]/estimator_array.mean[Full_R2_INDEX])**2+(trace_1.se/trace_1.mean)**2+(estimator_array.se[TR_RHO2_INDEX]/estimator_array.mean[TR_RHO2_INDEX])**2)
    return Stats(mean,se)

def stats_rdm_renyi_2(estimator_array):
    mean = -log(estimator_array.mean[RDM_R2_INDEX]/(estimator_array.mean[RDM_TR1_INDEX]*estimator_array.mean[RDM_TR2_INDEX]))/log(2)
    se = (1/log(2))*scipy.sqrt((estimator_array.se[RDM_R2_INDEX]/estimator_array.mean[RDM_R2_INDEX])**2+(estimator_array.se[RDM_TR1_INDEX]/estimator_array.mean[RDM_TR1_INDEX])**2+(estimator_array.se[RDM_TR2_INDEX]/estimator_array.mean[RDM_TR2_INDEX])**2)
    return Stats(mean,se)

def extract_data(data_files, estimtor_col_no, start_averaging):
    data = {}
    beta_loop_count = 1 
    # ignore lines beginning with # as the first non-whitespace character.
    comment = '^ *#'
    # start of data table
    start_regex = '^ # iterations'
    beta_loop_regex = '^ # Resetting beta'
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
                    if beta_loop_count >= start_averaging:
                        words = line.split()
                        if int(words[PARTICLES_COL]) == 0:
                             have_data = False
                        else:
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
                elif re.search(beta_loop_regex, line):
                    beta_loop_count = beta_loop_count + 1
                    have_data = True
            elif re.search(start_regex, line):
                have_data = True
            elif re.search(tau_regex, line):
                # get tau
                tau = float(re.search(tau_regex, line).group())

        f.close()

    return (data, tau)

def reweight_data(data, tau, l):
    '''Multiply data by necessary factors to remove population control bias.'''
    # Put the keys into a list.
    beta_values = list(data)
    # Sort the list into order, from smallest to largest.
    beta_values = sorted(beta_values)
    counter = 0
    max_index_numerators = len(data[0].numerators[0])
    max_index_loops = len(data[0].shift)
    
    for x in range(0,max_index_loops):
        accumulated_weight = 1.0
        counter = 0
        for beta in beta_values:
            # If more iterations have passed than the number of reweighting factors to be
            # used, then remove the oldest factors of e(-tau*shift) from accumulated_weight
            # so that only the last l factors of e(-tau*shift) will be included.
            if counter >= l:
                a = counter-l
                b = data[beta_values[a]].shift[x]
                accumulated_weight = accumulated_weight*exp(tau*b)
            b = data[beta].shift[x]
            accumulated_weight = accumulated_weight*exp(-tau*b)
            counter = counter+1
            data[beta].Tr_rho[x] = data[beta].Tr_rho[x]*accumulated_weight
            for y in range(0,max_index_numerators):
                data[beta].numerators[x][y] = data[beta].numerators[x][y]*accumulated_weight
    
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
        if Full_R2_is_present:
            stats[beta].full_renyi_2 = stats_full_renyi_2(stats[beta].estimators, stats[beta].Tr_rho)
        if RDM_R2_is_present:
            stats[beta].rdm_renyi_2 = stats_rdm_renyi_2(stats[beta].estimators)

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

def print_stats(stats, estimator_headings, trace=False,  shift=False, with_spline=False, with_heat_capacity=False):
    energies = []
#    energies_squared = []
    betas = []
    weights = []
#    weights_squared = []
    if with_spline:
        for beta in sorted(stats.iterkeys()):
            data = stats[beta]
            weights.append(1./data.estimators.se[TR_HRHO_INDEX])
#            weights_squared.append(1./data.estimators.se[TR_H2RHO_INDEX])
            betas.append(beta)
            energies.append(data.estimators.mean[TR_HRHO_INDEX]) 
#            energies_squared.append(data.estimators.mean[TR_H2RHO_INDEX]) 
        energy_fit, specific_heats = calculate_specific_heat(energies, betas, weights)
            
    print '#           beta    ',
    if trace:
        print 'trace     trace s.e.',
    if shift:
        print 'shift           shift s.e.    ',
    for i in range(0,len(estimator_headings)):
        if not (i == RDM_R2_INDEX or i == RDM_TR1_INDEX or i == RDM_TR2_INDEX or i == TR_RHO2_INDEX or i == Full_R2_INDEX):
            print estimator_headings[i]+'            s.e.    ',
    if H2_is_present and H_is_present and with_heat_capacity:
       print 'Stoch. Spec. Heat  s.e.    '
    if Full_R2_is_present:
       print ' Full Renyi 2    Full Renyi 2 s.e.     '
    if RDM_R2_is_present:
       print ' RDM Renyi 2     RDM Renyi 2 s.e.     '
    if with_spline:    
       print '  Energy Fit    Spline HC'
    print
    counter = 0
    for beta in sorted(stats.iterkeys()):

        data = stats[beta]
        print '%16.8f' % (beta) ,
        if trace:
            print '%16.8f%16.8f' % (data.Tr_rho.mean, data.Tr_rho.se) ,
        if shift:
            print '%16.8f%16.8f' % (data.shift.mean, data.shift.se) ,
        for i in range(0,len(data.estimators.mean)):
            if not (i == RDM_R2_INDEX or i == RDM_TR1_INDEX or i == RDM_TR2_INDEX or i == TR_RHO2_INDEX or i == Full_R2_INDEX):
                print '%16.8f%16.8f' % (data.estimators.mean[i], data.estimators.se[i]) ,
        if H2_is_present and H_is_present and with_heat_capacity:
            print '%16.8f%16.8f' % (data.stochastic_specific_heat.mean, data.stochastic_specific_heat.se) ,
        if Full_R2_is_present:
            print '%16.8f%16.8f' % (data.full_renyi_2.mean, data.full_renyi_2.se) ,
        if RDM_R2_is_present:
            print '%16.8f%16.8f' % (data.rdm_renyi_2.mean, data.rdm_renyi_2.se) ,
        if with_spline:
            print '%16.8f%16.8f' % (energy_fit[counter], specific_heats[counter]) ,
        counter = counter + 1
        print

def plot_stats(stats, data_file, trace=False, shift=False, spline=False, heat_capacity=False):

    if USE_MATPLOTLIB:
        save_name = data_file.split()
        save_temp = data_file.split(r'/') 
        save_name = save_temp[len(save_temp)-1].split('.',1)[0]

        fig_width_pt = 246.0  # Get this from LaTeX using \showthe\columnwidth
        inches_per_pt = 1.0/72.27               # Convert pt to inch
        golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
        fig_width = fig_width_pt*inches_per_pt  # width in inches
        fig_height = fig_width*golden_mean      # height in inches
        fig_size =  [fig_width,fig_height]

        pyplot.rcParams.update(
            {'legend.borderpad':0.3,
             'legend.labelspacing':0.3,
             'legend.fontsize':'medium',
             'text.usetex':True,
             'backend': 'ps',
             'axes.labelsize': 8,
             'text.fontsize': 8,
             'legend.fontsize': 8,
             'xtick.labelsize': 8,
             'ytick.labelsize': 8,
             'figure.figsize': fig_size}
        )

        beta = sorted(stats.keys())
        beta_max = max(beta)
        # find nice frequency to plot points
        N = int(len(beta)/25.)
        if shift:
            # data
            data = [stats[b].shift.mean for b in beta] 
            # errorbars
            err = [stats[b].shift.se for b in beta] 
            pyplot.gca().cla()
            pyplot.plot(beta, data, c='blue')
            pyplot.errorbar(beta[::N], data[::N], yerr=err[::N], linewidth = 0., c='blue')
            ax = pyplot.gca()
            ax.xaxis.set_major_formatter(pyplot.FormatStrFormatter(r'$%3i$'))
            ax.yaxis.set_major_formatter(pyplot.FormatStrFormatter(r'$%2.1f$'))
            ax.xaxis.set_major_locator(pyplot.MaxNLocator(5))
            ax.yaxis.set_major_locator(pyplot.MaxNLocator(4))
            ax.set_xlim((0,beta_max))
            y_min = min(data) - 0.05*(max(data) - min(data))
            y_max = max(data) + 0.05*(max(data) - min(data))
            ax.set_ylim((y_min,y_max))
            ax.set_xlabel(r'$\beta J$')
            ax.set_ylabel(r'$\textnormal{Shift}/J$')
            pyplot.subplots_adjust(left=0.15, bottom=0.17, right=0.95)
            pyplot.savefig(save_name+'_shift.eps', dpi=(300))
            pyplot.savefig(save_name+'_shift.pdf')
            
        if trace:

            # data
            data = [stats[b].Tr_rho.mean/1000000. for b in beta]
            # errorbars
            err = [stats[b].Tr_rho.se/1000000. for b in beta]
            pyplot.gca().cla()
            pyplot.plot(beta, data, c='blue')
            pyplot.errorbar(beta[::N], data[::N], yerr=err[::N], linewidth = 0., c='blue')            
            ax = pyplot.gca()
            ax.xaxis.set_major_formatter(pyplot.FormatStrFormatter(r'$%3i$'))
            ax.yaxis.set_major_formatter(pyplot.FormatStrFormatter(r'$%2.1f$'))
            ax.xaxis.set_major_locator(pyplot.MaxNLocator(5))
            ax.yaxis.set_major_locator(pyplot.MaxNLocator(4))            
            ax.set_xlim((0,beta_max))
            y_min = min(data) - 0.05*(max(data) - min(data))
            y_max = max(data) + 0.05*(max(data) - min(data))
            ax.set_ylim((y_min,y_max))            
            ax.set_xlabel(r'$\beta J$')
            ax.set_ylabel(r'$\sum_{i}\rho_{ii}/1\times 10^6 $')
            pyplot.subplots_adjust(left=0.15, bottom=0.17, right=0.95)
            pyplot.savefig(save_name+'_trace.eps', dpi=(300))
            pyplot.savefig(save_name+'_trace.pdf')

        if H2_is_present and H_is_present and  heat_capacity:
      
            # data
            data = [stats[b].stochastic_specific_heat.mean for b in beta]
            # errorbars
            err = [stats[b].stochastic_specific_heat.se for b in beta]
            pyplot.gca().cla()
            pyplot.plot(beta, data, c='blue')
            pyplot.errorbar(beta[::N], data[::N], yerr=err[::N], linewidth = 0., c='blue')
            ax = pyplot.gca()

            ax.set_xlim((0,beta_max))
            y_min = min(data) - 0.05*(max(data) - min(data))
            y_max = max(data) + 0.05*(max(data) - min(data))
            ax.set_ylim((y_min,y_max))            
            ax.set_xlabel(r'$\beta J$')
            ax.set_ylabel(r'$\langle C_{h} \rangle /k_B$')
            pyplot.subplots_adjust(left=0.15, bottom=0.17, right=0.95)
            pyplot.savefig(save_name+'_stochastic_specific_heat.eps', dpi=(300))
            pyplot.savefig(save_name+'_stochastic_specific_heat.pdf')

        if H2_is_present:
            
            # data
            data = [stats[b].estimators.mean[TR_H2RHO_INDEX] for b in beta]
            # errorbars
            err = [stats[b].estimators.se[TR_H2RHO_INDEX] for b in beta]
            pyplot.gca().cla()
            pyplot.plot(beta, data, c='blue')
            pyplot.errorbar(beta[::N], data[::N], yerr=err[::N], linewidth = 0., c='blue')
            ax = pyplot.gca()
            ax.xaxis.set_major_formatter(pyplot.FormatStrFormatter(r'$%3i$'))
            ax.yaxis.set_major_formatter(pyplot.FormatStrFormatter(r'$%2.1f$'))
            ax.xaxis.set_major_locator(pyplot.MaxNLocator(5))
            ax.yaxis.set_major_locator(pyplot.MaxNLocator(4))
            ax.set_xlim((0,beta_max))
            y_min = min(data) - 0.05*(max(data) - min(data))
            y_max = max(data) + 0.05*(max(data) - min(data))
            ax.set_ylim((y_min,y_max)) 
            ax.set_xlabel(r'$\beta J$')
            ax.set_ylabel(r'$\langle H^2 \rangle /J^2$')
            pyplot.subplots_adjust(left=0.15, bottom=0.17, right=0.95)
            pyplot.savefig(save_name+'_energy_squared.eps', dpi=(300))
            pyplot.savefig(save_name+'_energy_squared.pdf')

        if H_is_present:
         
            # data
            data = [stats[b].estimators.mean[TR_HRHO_INDEX] for b in beta]
            # errorbars
            err = [stats[b].estimators.se[TR_HRHO_INDEX] for b in beta]
            pyplot.gca().cla()
            pyplot.plot(beta, data, c='blue')
            pyplot.errorbar(beta[::N], data[::N], yerr=err[::N], linewidth = 0., c='blue')
            ax = pyplot.gca()
            ax.xaxis.set_major_formatter(pyplot.FormatStrFormatter(r'$%3i$'))
            ax.yaxis.set_major_formatter(pyplot.FormatStrFormatter(r'$%2.1f$'))
            ax.xaxis.set_major_locator(pyplot.MaxNLocator(5))
            ax.yaxis.set_major_locator(pyplot.MaxNLocator(4))
            ax.set_xlim((0,beta_max))
            y_min = min(data) - 0.05*(max(data) - min(data))
            y_max = max(data) + 0.05*(max(data) - min(data))
            ax.set_ylim((y_min,y_max))            
            ax.set_xlabel(r'$\beta J$')
            ax.set_ylabel(r'$\langle H \rangle /J$')
            pyplot.subplots_adjust(left=0.15, bottom=0.17, right=0.95)
            pyplot.savefig(save_name+'_energy.eps', dpi=(300))
            pyplot.savefig(save_name+'_energy.pdf')
  
        if S_is_present: 
            # data
            data = [stats[b].estimators.mean[TR_SRHO_INDEX] for b in beta]
            # errorbars
            err = [stats[b].estimators.se[TR_SRHO_INDEX] for b in beta]
            pyplot.gca().cla()
            pyplot.plot(beta, data, c='blue')
            pyplot.errorbar(beta[::N], data[::N], yerr=err[::N], linewidth = 0., c='blue')            
            ax = pyplot.gca()
            ax.xaxis.set_major_formatter(pyplot.FormatStrFormatter(r'$%3i$'))
            ax.yaxis.set_major_formatter(pyplot.FormatStrFormatter(r'$%2.1f$'))
            ax.xaxis.set_major_locator(pyplot.MaxNLocator(5))
            ax.yaxis.set_major_locator(pyplot.MaxNLocator(4))
            ax.set_xlim((0,beta_max))
            y_min = min(data) - 0.05*(max(data) - min(data))
            y_max = max(data) + 0.05*(max(data) - min(data))
            ax.set_ylim((y_min,y_max))            
            ax.set_xlabel(r'$\beta J$')
            ax.set_ylabel(r'$\langle M /rangle$')
            pyplot.subplots_adjust(left=0.15, bottom=0.17, right=0.95)
            pyplot.savefig(save_name+'_staggerred_magnetisation.eps', dpi=(300))
            pyplot.savefig(save_name+'_staggered_magnetisation.pdf')

        if C_is_present: 
            # data
            data = [stats[b].estimators.mean[TR_CRHO_INDEX] for b in beta]
            # errorbars
            err = [stats[b].estimators.se[TR_CRHO_INDEX] for b in beta]
            pyplot.gca().cla()
            pyplot.plot(beta, data, c='blue')
            pyplot.errorbar(beta[::N], data[::N], yerr=err[::N], linewidth = 0., c='blue')
            ax = pyplot.gca()
            ax.xaxis.set_major_formatter(pyplot.FormatStrFormatter(r'$%3i$'))
            ax.yaxis.set_major_formatter(pyplot.FormatStrFormatter(r'$%2.1f$'))
            ax.xaxis.set_major_locator(pyplot.MaxNLocator(5))
            ax.yaxis.set_major_locator(pyplot.MaxNLocator(4))
            ax.set_xlim((0,beta_max))
            y_min = min(data) - 0.05*(max(data) - min(data))
            y_max = max(data) + 0.05*(max(data) - min(data))
            ax.set_ylim((y_min,y_max))            
            ax.set_xlabel(r'$\beta J$')
            ax.set_ylabel(r'$\langle C \rangle$')
            pyplot.subplots_adjust(left=0.15, bottom=0.17, right=0.95)
            pyplot.savefig(save_name+'_spin_correlation.eps', dpi=(300))
            pyplot.savefig(save_name+'_spin_correlation.pdf')

        if M2_is_present:

            # data
            data = [stats[b].estimators.mean[TR_M2RHO_INDEX] for b in beta]
            # errorbars
            err = [stats[b].estimators.se[TR_M2RHO_INDEX] for b in beta]
            pyplot.gca().cla()
            pyplot.plot(beta, data, c='blue')
            pyplot.errorbar(beta[::N], data[::N], yerr=err[::N], linewidth = 0., c='blue')
            ax = pyplot.gca()
            ax.xaxis.set_major_formatter(pyplot.FormatStrFormatter(r'$%3i$'))
            ax.yaxis.set_major_formatter(pyplot.FormatStrFormatter(r'$%2.1f$'))
            ax.xaxis.set_major_locator(pyplot.MaxNLocator(5))
            ax.yaxis.set_major_locator(pyplot.MaxNLocator(4))
            ax.set_xlim((0,beta_max))
            y_min = min(data) - 0.05*(max(data) - min(data))
            y_max = max(data) + 0.05*(max(data) - min(data))
            ax.set_ylim((y_min,y_max))            
            ax.set_xlabel(r'$\beta J$')
            ax.set_ylabel(r'$\langle M^2 \rangle$')
            pyplot.subplots_adjust(left=0.15, bottom=0.17, right=0.95)
            pyplot.savefig(save_name+'_magnetisation_squared.eps', dpi=(300))
            pyplot.savefig(save_name+'_magnetisation_squared.pdf')

 
def parse_options(args):

    parser = optparse.OptionParser(usage = __doc__)
    parser.add_option('--plot', action='store_true', default=False, help='Plot averaged data.')
    parser.add_option('--with-shift', action='store_true', dest='with_shift', default=False, help='Analyse shift data.')
    parser.add_option('--without-shift', action='store_false', dest='with_shift', help='Do not analyse shift data.  Default.')
    parser.add_option('--with-trace', action='store_true', dest='with_trace', default=False, help='Analyse trace data.  Default.')
    parser.add_option('--without-trace', action='store_false', dest='with_trace', help='Do not analyse trace data.')
    parser.add_option('--with-spline', action='store_true', dest='with_spline', default=False, help='Calculate heat capacity using a spline fit of H2 and H')
    parser.add_option('--without-spline', action='store_false', dest='with_spline', help='Do not calcualate heat capacity using spline fit. Default')
    parser.add_option('--with-heat-capacity', action='store_true', dest='with_heat_capacity', default=False, help='Calculate the stochastic heat capacity.')
    parser.add_option('--without-heat-capacity', action='store_false', dest='with_heat_capacity', help='Do not calculate teh stochastic heat capacity')
    parser.add_option('-b','--remove_bias', action='store_true', dest='remove_bias', default=False, help='Undo the biasing effect of the shift by reweighting the data.')
    parser.add_option('-a','--averaging_length', dest='averaging_length', type='int', default=100, help='When removing the bias by reweighting, this gives the number of previous iterations over which the weights are accumulated. Default: %default.')
    parser.add_option("-s", "--start-averaging", action="store", dest="start_averaging", type="int")
    (options, filenames) = parser.parse_args(args)
    
    if len(filenames) == 0:
        parser.print_help()
        sys.exit(1)

    return (options, filenames)

if __name__ == '__main__':
    (options, data_files) = parse_options(sys.argv[1:])
    estimator_headings, estimtor_col_no = get_estimator_headings(data_files)
    (data,tau) = extract_data(data_files, estimtor_col_no, options.start_averaging)
    if options.remove_bias:
        data = reweight_data(data,tau,options.averaging_length)
    stats = get_data_stats(data)
    print_stats(stats, estimator_headings, options.with_trace, options.with_shift, options.with_spline, options.with_heat_capacity)
    if options.plot:
        plot_stats(stats, data_files[0], options.with_trace, options.with_shift, options.with_spline, options.with_heat_capacity)
