'''pyhande provides powerful abstractions for analysing HANDE calculations.

HANDE includes many scripts for common analysis tasks which are (typically) thin
wrappers around pyhande.  More complicated data analysis or examining large
numbers of output files can be easily performed by using pyhande directly from a
python interpreter or a custom script.

pyhande can extract data from output produced by HANDE and perform a variety of
data analysis tasks on the data obtained.  See the documentation for each 
submodule for more details.
'''

# For convenience, import all submodules so the user need only import pyhande.
import pyhande.analysis
import pyhande.canonical
import pyhande.extract
import pyhande.lazy
import pyhande.utils
import pyhande.weight
