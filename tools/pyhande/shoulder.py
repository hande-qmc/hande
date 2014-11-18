''' Estimate the shoulder (plateau) from a (CCMC) FCIQMC '''

import pandas as pd
import numpy as np

def shoulder_estimator(data, total_key='# H psips', ref_key='N_0'):
    ''' Estimate the shoulder (plateau) from a (CCMC) FCIQMC calulation.
The shoulder esitmator is deffined to be the ten points with the smallest
proportion of the poulation on the reference (excluding points when the 
population drops bellow 10 excips (psips). The shoulder height is the total
population at this point.

Parameters
----------
data : :class:`pandas.DataFrame`
    HANDE QMC data. The function pyhande.extract.extract_data_sets can 
    be used to extract this from a HANDE outputfile.
total_key : string
    column name in reblock_data containing the total number of psips.
ref_key : string
    column name in reblock_data containing the number of psips on the
    reference determinant.

Returns
-------
platue_data : :class:`pandas.DataFrame`
    An estimate of the shoulder (platue) from a FCIQMC (CCMC) calulation,
    along with the associated standerd error.
'''
    platue_data = []

    #Set data less than 10 to 10.
    #We explicity have to copy the data frame otherwise we will edit the actual data frame
    data_1 = data.copy()
    data_1[ref_key][data_1[ref_key]<10] = np.ones(len(data_1[ref_key][data_1[ref_key]<10]))*10
    data_1['Shoulder'] = data_1[total_key]/data_1[ref_key]
    sorted_data = data_1.sort('Shoulder')
    platue_data.append([sorted_data[-10:]['Shoulder'].mean(), sorted_data[-10:]['Shoulder'].sem()])
    platue_data.append([sorted_data[-10:][total_key].mean(), sorted_data[-10:][total_key].sem()])
    platue_data = pd.DataFrame(data=platue_data, columns=['mean', 'standard error'],
                               index=['shoulder estimator', 'shoulder height'])
    return platue_data
