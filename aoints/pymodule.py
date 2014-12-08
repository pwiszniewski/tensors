import psi4
import re
import os
import inputparser
import math
import warnings
from driver import *
from wrappers import *
from molutil import *
import p4util
from psiexceptions import *


def run_aoints(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    aoints can be called via :py:func:`~driver.energy`. For post-scf plugins.

    >>> energy('aoints')

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    # Your plugin's psi4 run sequence goes here
    psi4.set_global_option('BASIS', 'sto-3g')
    psi4.set_local_option('AOINTS', 'PRINT', 1)
    scf_helper(name, **kwargs)
    returnvalue = psi4.plugin('aoints.so')
    psi4.set_variable('CURRENT ENERGY', returnvalue)


# Integration with driver routines
procedures['energy']['aoints'] = run_aoints


def exampleFN():
    # Your Python code goes here
    pass
