.. _state_histogram_table:

state_histogram options
=======================

The ``state_histogram`` table contains options used to control the options for
calculating the state histograms during a :ref:`fciqmc` or :ref:`dmqmc`
simulation.

This is a generalization of the state histograms original proposed by
Cleland and coworkers for FCIQMC: [Cleland12]_. The generalization is to treat
DMQMC in addition to FCIQMC. In the case of FCIQMC the behavior is nearly
identical to that of the original algorithm.

In FCIQMC, the histograms are calculate based on the excitation levels relative
to the reference determinant. In DMQMC, we add an additional excitation index
based on the excitation level between the two determinants used to label the
density matrix site.

``report_frequency``
    type: integer.

    Optional. Default: :ref:`nreports<qmc_table>`

    The frequency in report cycles to calculate and report the state histograms.
    A histogram is always reported at the beginning and end of a calculation
    regardless of whether it falls within the frequency provided.

``nbins``
    type: integer.

    Optional. Default: 5

    The number of histogram bins to use per decade of walker population.
    As an example, if we have a walker population 10 to 100, there are bins
    from: :math:`\left[10, \sim 15.85\right)`, :math:`\left[15.85, \sim 25.12\right)`,
    :math:`\left[25.12, \sim 39.81\right)`, :math:`\left[39.81, \sim 63.10\right)`, and
    :math:`\left[63.10, \sim 100\right)`.

``decades``
    type: integer.

    Optional. Default: 3

    The number of decades of particle populations to include in the histogram
    bins past the decades of :math:`\lfloor \log_{10}(N_w) \rfloor`. Where
    :math:`N_w` is the simulations ``target_population`` set in the :ref:`qmc_table`.

    .. note::

        If this value is too small, and the simulations population grows well
        pass the target population, the histogram population bins may not cover the number
        of particle on a single site. In this case the simulation will terminate
        with an error message to prevent an out of bounds index.

``skip_memory_check``
    type: boolean.

    Optional. Default: false.

    Controls the memory check performed when collecting state histograms.
    Normally an estimate is made for the memory cost to store all the state
    histogram files, if the estimate exceeds 1 gigabyte the calculation will
    halt. Setting this flag to true will ignore the memory check.
