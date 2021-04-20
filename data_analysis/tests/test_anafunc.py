import unittest
import time
import sys
sys.path.append("C:\\Users\\e.zippo\\Desktop\\Università\\Tesi\\codice_lattice_su2\\lattice-su2")

import numpy as np
from matplotlib import pyplot as pyp
from scipy.optimize import curve_fit

import data_analysis.analysis_functions as ana_func 

class TestAnaFunc(unittest.TestCase):
    ''' test of analisys_functions.py '''
    def test_blocking(self):
        sample = np.random.normal(0,1,1000)
        block_new = ana_func.sigma_blocking(sample, 3)
        block_old = ana_func.sigma_blocking_old(sample, 3)
        self.assertAlmostEqual(block_new,block_old)
    
if __name__=='__main__':
    unittest.main()
    
