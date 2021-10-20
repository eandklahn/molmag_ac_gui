# try to import known packages
# this is to check whether all needed packages are there
# before anything is invoked, so that the GUI does not
# raise an error unexpectedly.
import numpy
import scipy
import pandas
import names
import matplotlib
import lmfit
import PyQt5
import PIL

# own packages
from . import ac_gui
from . import process_ac

__all__ = [
           'process_ac',
           'ac_gui'
          ]

__version__ = '0.2.6'
__author__ = 'Emil Andreasen Klahn'
__author_email__ = 'eklahn@chem.au.dk'