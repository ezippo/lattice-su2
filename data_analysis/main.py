import sys
sys.path.append("C:\\Users\\e.zippo\\Desktop\\Università\\Tesi\\codice_lattice_su2\\lattice-su2")

import numpy as np
from matplotlib import pyplot as pyp
from scipy.optimize import curve_fit

from data_analysis import analysis_functions as ana_func

if __name__=='__main__':
    dati = np.loadtxt("C:\\Users\\e.zippo\\Desktop\\Università\\Tesi\\codice_lattice_su2\\su2-1.0.0\\dati.dat")
    n_obs = 3
    mu1 = ana_func.average(dati[:,n_obs])
    sigma1 = ana_func.sigma_bootstrap(dati[:,n_obs], 3, 300)
    print(mu1, sigma1)
    
    ana_func.plot_sampletrend(dati[:,n_obs])
    pyp.show()