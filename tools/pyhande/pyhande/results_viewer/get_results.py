"""Helper functions to run analysis and get results object."""
from typing import List, Tuple, Union
from pyhande.data_preparing.hande_ccmc_fciqmc import PrepHandeCcmcFciqmc
from pyhande.extracting.extractor import Extractor
from pyhande.error_analysing.blocker import Blocker
from pyhande.error_analysing.hybrid_ana import HybridAna
from pyhande.helpers.simple_callables import RaiseValueError
from pyhande.results_viewer.results import Results
from pyhande.results_viewer.results_ccmc_fciqmc import ResultsCcmcFciqmc


def define_objects_common(merge_type: str = 'uuid', analyser: str = 'blocking',
                          start_its: Union[List[int], str] = 'blocking'
                          ) -> Tuple[Extractor, PrepHandeCcmcFciqmc,
                                     Union[Blocker, HybridAna]]:
    """Create extractor, preparator and analyser with common options.

    Parameters
    ----------
    merge_type : str, optional
        how to do merge, 'uuid', 'legacy' or 'no. Note that this is
        different to fuller options when instantiating extractor object
        directly, by default 'uuid'.
    analyser : str, optional
        'blocking' for doing reblocking or 'hybrid',
        by default 'blocking'
    start_its : Union[List[int], str], optional
        Either list of integer for start iterations or 'blocking' or
        'hybrid', defining find starting iteration function to use.
        by default 'blocking'

    Returns
    -------
    Tuple[Extractor, PrepHandeCcmcFciqmc, Union[Blocker, HybridAna]]
        Instantiated objects for extracting, preparing and analysing
        data.
    """
    extra = Extractor(merge={'type': merge_type})
    prep = PrepHandeCcmcFciqmc()
    select_ana = {'blocking': Blocker.inst_hande_ccmc_fciqmc,
                  'hybrid': HybridAna.inst_hande_ccmc_fciqmc}
    ana = select_ana.get(analyser, RaiseValueError("Choose valid analyser."))(
        start_its=start_its)
    return (extra, prep, ana)


def analyse_data(out_files: List[str], extractor: Extractor,
                 preparator: PrepHandeCcmcFciqmc = None,
                 analyser: Union[Blocker, HybridAna] = None
                 ) -> Union[Results, ResultsCcmcFciqmc]:
    """Execute objects to extract data, prepare and analyse it.

    Parameters
    ----------
    out_files: List[str]
        Output files with data to extract, prepare and analyse.
    extractor: Extractor
        Instance to extract data from files.
    preparator: PrepHandeCcmcFciqmc
        Instance to prepare data, e.g. calculate inst. proj. energy or
        deal with complex/replica tricks. The default is None.
    analyser: Union[Blocker, HybridAna]
        Instance to analyse data, e.g. blocking. The default is None.

    Returns
    -------
    Union[Results, ResultsCcmcFciqmc]
        Results object to view and further analyse results.
    """
    extractor.exe(out_files)
    if extractor.all_ccmc_fciqmc:
        if preparator:
            preparator.exe(extractor.data)
            if analyser:
                analyser.exe(preparator.data, preparator.observables)
        return ResultsCcmcFciqmc(extractor, preparator, analyser)
    return Results(extractor)


def get_results(out_files: List[str], merge_type: str = 'uuid',
                analyser: str = 'blocking',
                start_its: Union[List[int], str] = 'blocking'
                ) -> Union[Results, ResultsCcmcFciqmc]:
    """Lazy function to combine defining objects and executing them.

    Parameters
    ----------
    see define_objects_common and analyse_data

    Returns
    -------
    see analyse_data
    """
    return analyse_data(
        out_files, *define_objects_common(merge_type, analyser, start_its))
