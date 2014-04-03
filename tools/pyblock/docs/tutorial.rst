
pyblock tutorial
================

The estimate of the standard error of a set of data assumes that the
data points are completely independent. If this is not true, then
naively calculating the standard error of the entire data set can give a
substantial underestimate of the true error. This arises in, for
example, Monte Carlo simulations where the state at one step depends
upon the state at the previous step. Data calculated from the stochastic
state hence has **serial correlations**.

A simple way to remove these correlations is to repeatedly average
neighbouring pairs of data points and calculate the standard error on
the new data set. As no data is discarded in this process (assuming the
data set contains :math:`2^n` values), the error estimate should remain
approximately constant if the data is truly independent.

:mod:``pyblock`` is a python module for performing this reblocking
analysis.

Normally correlated data comes from an experiment or simulation but
we'll use randomly generated data which is serially correlated in order
to show how :mod:\`pyblock\` works.

.. code:: python

    # Running in an ipython notebook (hence some semicolons to suppress unwanted output below).
    %pylab inline

.. parsed-literal::

    Populating the interactive namespace from numpy and matplotlib


.. code:: python

    import numpy
    def corr_data(N, L):
        '''Generate random correlated data containing 2^N data points.  
        Randon data is convolved over a 2^L/10 length to give the correlated signal.'''
        return numpy.convolve(numpy.random.randn(2**N), numpy.ones(2**L)/10, 'same')
    rand_data = corr_data(16, 6)
.. code:: python

    plot(rand_data);


.. image:: tutorial_files/tutorial_4_0.png


If we zoom in, we can clearly see that neighbouring data points do not
immediately appear to be independent:

.. code:: python

    plot(rand_data[:1000]);
    plot(rand_data[40000:41000]);


.. image:: tutorial_files/tutorial_6_0.png


``pyblock`` can perform a reblocking analysis to get a better estimate
of the standard error of the data set:

.. code:: python

    import pyblock
    reblock_data = pyblock.blocking.reblock(rand_data)
    for reblock_iter in reblock_data:
        print(reblock_iter)

.. parsed-literal::

    BlockTuple(block=0, ndata=65536, mean=array(0.02796172405542975), cov=array(0.6006764511115371), std_err=array(0.0030274734123228956), std_err_err=array(8.36235668947543e-06))
    BlockTuple(block=1, ndata=32768, mean=array(0.02796172405542985), cov=array(0.5956967293152431), std_err=array(0.0042637098254553986), std_err_err=array(1.6655370648411603e-05))
    BlockTuple(block=2, ndata=16384, mean=array(0.02796172405542986), cov=array(0.5881998793984137), std_err=array(0.005991733600996738), std_err_err=array(3.310097471080922e-05))
    BlockTuple(block=3, ndata=8192, mean=array(0.027961724055429855), cov=array(0.5745686682483024), std_err=array(0.008374829961603943), std_err_err=array(6.543235287079677e-05))
    BlockTuple(block=4, ndata=4096, mean=array(0.027961724055429925), cov=array(0.5480484251400315), std_err=array(0.011567233249310442), std_err_err=array(0.00012781668279599132))
    BlockTuple(block=5, ndata=2048, mean=array(0.027961724055429883), cov=array(0.4936841949486424), std_err=array(0.015526001926277311), std_err_err=array(0.00024265302879358065))
    BlockTuple(block=6, ndata=1024, mean=array(0.02796172405542991), cov=array(0.40100754763760316), std_err=array(0.019789111481818653), std_err_err=array(0.00043749538930113217))
    BlockTuple(block=7, ndata=512, mean=array(0.027961724055429887), cov=array(0.24257748878302304), std_err=array(0.02176658351187301), std_err_err=array(0.0006808709727874168))
    BlockTuple(block=8, ndata=256, mean=array(0.02796172405542989), cov=array(0.1417545351603556), std_err=array(0.02353143967057985), std_err_err=array(0.0010419896625227117))
    BlockTuple(block=9, ndata=128, mean=array(0.02796172405542988), cov=array(0.0768776702819086), std_err=array(0.02450728053206661), std_err_err=array(0.0015377235437994572))
    BlockTuple(block=10, ndata=64, mean=array(0.027961724055429887), cov=array(0.033939374250020415), std_err=array(0.023028302643846095), std_err_err=array(0.002051524254576491))
    BlockTuple(block=11, ndata=32, mean=array(0.027961724055429897), cov=array(0.014287585058883745), std_err=array(0.02113023977833941), std_err_err=array(0.0026835431353935822))
    BlockTuple(block=12, ndata=16, mean=array(0.027961724055429887), cov=array(0.006046924816731885), std_err=array(0.01944049384778439), std_err_err=array(0.003549332336490638))
    BlockTuple(block=13, ndata=8, mean=array(0.02796172405542989), cov=array(0.0019730090976917227), std_err=array(0.015704334981509575), std_err_err=array(0.0041971600705669795))
    BlockTuple(block=14, ndata=4, mean=array(0.027961724055429887), cov=array(0.0011946223767190995), std_err=array(0.017281654844943956), std_err_err=array(0.007055206046834906))
    BlockTuple(block=15, ndata=2, mean=array(0.027961724055429887), cov=array(0.0005757022817790697), std_err=array(0.016966176378003822), std_err_err=array(0.011996898367693519))


The standard error of the original data set is clearly around 8 times
too small. Note that the standard error of the last few reblock
iterations fluctuates substantially---this is simply because of the
small number of data points at those iterations.

In addition to the mean and standard error at each iteration, the
covariance and an estimate of the error in the standard error are also
calculated. Each tuple also contains the number of data points used at
the given reblock iteration.

``pyblock`` can also suggest the reblock iteration at which the standard
error has converged (i.e. the iteration at which the serial correlation
has been removed and every data point is truly independent).

.. code:: python

    opt = pyblock.blocking.find_optimal_block(len(rand_data), reblock_data)
    print(opt)
    print(reblock_data[opt[0]])

.. parsed-literal::

    [10]
    BlockTuple(block=10, ndata=64, mean=array(0.027961724055429887), cov=array(0.033939374250020415), std_err=array(0.023028302643846095), std_err_err=array(0.002051524254576491))


Whilst the above uses just a single data set, ``pyblock`` is designed to
work on multiple data sets at once (e.g. multiple outputs from the same
simulation). In that case, different optimal reblock iterations might be
found for each data set. The only assumption is that the original data
sets are of the same length.

pandas integration
------------------

The core ``pyblock`` functionality is built upon ``numpy``. However, it
is more convenient to use the ``pandas``-based wrapper around
``pyblock.blocking``, not least because it makes working with multiple
data sets more pleasant.

.. code:: python

    import pandas as pd
    rand_data = pd.Series(rand_data)
.. code:: python

    rand_data.head()



.. parsed-literal::

    0   -0.821195
    1   -0.891607
    2   -0.953600
    3   -1.030329
    4   -0.954974
    dtype: float64



.. code:: python

    (data_length, reblock_data, covariance) = pyblock.pd_utils.reblock(rand_data)
.. code:: python

    # number of data points at each reblock iteration
    data_length



.. parsed-literal::

    reblock
    0          65536
    1          32768
    2          16384
    3           8192
    4           4096
    5           2048
    6           1024
    7            512
    8            256
    9            128
    10            64
    11            32
    12            16
    13             8
    14             4
    15             2
    Name: data length, dtype: int64



.. code:: python

    # mean, standard error and estimate of the error in the standard error at each reblock iteration
    # Note the suggested reblock iteration is already indicated.
    # pyblock names the data series 'data' if no name is provided in the pandas.Series/pandas.DataFrame.
    reblock_data



.. raw:: html

    <div style="max-height:1000px;max-width:1500px;overflow:auto;">
    <table border="1" class="dataframe">
      <thead>
        <tr>
          <th></th>
          <th colspan="4" halign="left">data</th>
        </tr>
        <tr>
          <th></th>
          <th>mean</th>
          <th>standard error</th>
          <th>standard error error</th>
          <th>optimal block</th>
        </tr>
        <tr>
          <th>reblock</th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0 </th>
          <td> 0.027962</td>
          <td> 0.003027</td>
          <td> 0.000008</td>
          <td>         </td>
        </tr>
        <tr>
          <th>1 </th>
          <td> 0.027962</td>
          <td> 0.004264</td>
          <td> 0.000017</td>
          <td>         </td>
        </tr>
        <tr>
          <th>2 </th>
          <td> 0.027962</td>
          <td> 0.005992</td>
          <td> 0.000033</td>
          <td>         </td>
        </tr>
        <tr>
          <th>3 </th>
          <td> 0.027962</td>
          <td> 0.008375</td>
          <td> 0.000065</td>
          <td>         </td>
        </tr>
        <tr>
          <th>4 </th>
          <td> 0.027962</td>
          <td> 0.011567</td>
          <td> 0.000128</td>
          <td>         </td>
        </tr>
        <tr>
          <th>5 </th>
          <td> 0.027962</td>
          <td> 0.015526</td>
          <td> 0.000243</td>
          <td>         </td>
        </tr>
        <tr>
          <th>6 </th>
          <td> 0.027962</td>
          <td> 0.019789</td>
          <td> 0.000437</td>
          <td>         </td>
        </tr>
        <tr>
          <th>7 </th>
          <td> 0.027962</td>
          <td> 0.021767</td>
          <td> 0.000681</td>
          <td>         </td>
        </tr>
        <tr>
          <th>8 </th>
          <td> 0.027962</td>
          <td> 0.023531</td>
          <td> 0.001042</td>
          <td>         </td>
        </tr>
        <tr>
          <th>9 </th>
          <td> 0.027962</td>
          <td> 0.024507</td>
          <td> 0.001538</td>
          <td>         </td>
        </tr>
        <tr>
          <th>10</th>
          <td> 0.027962</td>
          <td> 0.023028</td>
          <td> 0.002052</td>
          <td> &lt;---    </td>
        </tr>
        <tr>
          <th>11</th>
          <td> 0.027962</td>
          <td> 0.021130</td>
          <td> 0.002684</td>
          <td>         </td>
        </tr>
        <tr>
          <th>12</th>
          <td> 0.027962</td>
          <td> 0.019440</td>
          <td> 0.003549</td>
          <td>         </td>
        </tr>
        <tr>
          <th>13</th>
          <td> 0.027962</td>
          <td> 0.015704</td>
          <td> 0.004197</td>
          <td>         </td>
        </tr>
        <tr>
          <th>14</th>
          <td> 0.027962</td>
          <td> 0.017282</td>
          <td> 0.007055</td>
          <td>         </td>
        </tr>
        <tr>
          <th>15</th>
          <td> 0.027962</td>
          <td> 0.016966</td>
          <td> 0.011997</td>
          <td>         </td>
        </tr>
      </tbody>
    </table>
    <p>16 rows × 4 columns</p>
    </div>



.. code:: python

    # Covariance matrix is not so relevant for a single data set.
    covariance



.. raw:: html

    <div style="max-height:1000px;max-width:1500px;overflow:auto;">
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th></th>
          <th>data</th>
        </tr>
        <tr>
          <th>reblock</th>
          <th></th>
          <th></th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0 </th>
          <th>data</th>
          <td> 0.600676</td>
        </tr>
        <tr>
          <th>1 </th>
          <th>data</th>
          <td> 0.595697</td>
        </tr>
        <tr>
          <th>2 </th>
          <th>data</th>
          <td> 0.588200</td>
        </tr>
        <tr>
          <th>3 </th>
          <th>data</th>
          <td> 0.574569</td>
        </tr>
        <tr>
          <th>4 </th>
          <th>data</th>
          <td> 0.548048</td>
        </tr>
        <tr>
          <th>5 </th>
          <th>data</th>
          <td> 0.493684</td>
        </tr>
        <tr>
          <th>6 </th>
          <th>data</th>
          <td> 0.401008</td>
        </tr>
        <tr>
          <th>7 </th>
          <th>data</th>
          <td> 0.242577</td>
        </tr>
        <tr>
          <th>8 </th>
          <th>data</th>
          <td> 0.141755</td>
        </tr>
        <tr>
          <th>9 </th>
          <th>data</th>
          <td> 0.076878</td>
        </tr>
        <tr>
          <th>10</th>
          <th>data</th>
          <td> 0.033939</td>
        </tr>
        <tr>
          <th>11</th>
          <th>data</th>
          <td> 0.014288</td>
        </tr>
        <tr>
          <th>12</th>
          <th>data</th>
          <td> 0.006047</td>
        </tr>
        <tr>
          <th>13</th>
          <th>data</th>
          <td> 0.001973</td>
        </tr>
        <tr>
          <th>14</th>
          <th>data</th>
          <td> 0.001195</td>
        </tr>
        <tr>
          <th>15</th>
          <th>data</th>
          <td> 0.000576</td>
        </tr>
      </tbody>
    </table>
    <p>16 rows × 1 columns</p>
    </div>



We can also plot the convergence of the standard error estimate and
obtain a summary of the suggested data to quote:

.. code:: python

    pyblock.pd_utils.plot_reblocking(reblock_data);


.. image:: tutorial_files/tutorial_21_0.png


The standard error clearly converges to ~0.022. The suggested reblock
iteration (which uses a slightly conservative formula) is indicated by
the arrow on the plot.

.. code:: python

    pyblock.pd_utils.reblock_summary(reblock_data)



.. raw:: html

    <div style="max-height:1000px;max-width:1500px;overflow:auto;">
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>mean</th>
          <th>standard error</th>
          <th>standard error error</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>data</th>
          <td> 0.02796172</td>
          <td> 0.0230283</td>
          <td> 0.002051524</td>
        </tr>
      </tbody>
    </table>
    <p>1 rows × 3 columns</p>
    </div>



``pyblock`` also contains simple error propogation functions for
combining multiple noisy data sets and can handle multiple data sets at
once (contained either within a ``numpy`` array using
``pyblock.blocking`` or within a ``pandas.DataFrame``.
