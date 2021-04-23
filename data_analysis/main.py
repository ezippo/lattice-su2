import sys
sys.path.append("C:\\Users\\e.zippo\\Desktop\\Università\\Tesi\\codice_lattice_su2\\lattice-su2\\data_analysis")

import numpy as np
from matplotlib import pyplot as pyp
from scipy.optimize import curve_fit

from data_analysis import analysis_functions as ana_func

def plaquette_tot(dati):
    plaq = np.array([ dati[:,0], dati[:,1] ])
    plaq = np.reshape(plaq, len(plaq[0])*2)
    return plaq

def plaquette_adj_tot(dati):
    if len(dati)<5:
        return 0
    plaq = np.array([ dati[:,2], dati[:,3] ])
    plaq = np.reshape(plaq, len(plaq[0])*2)
    return plaq
    
def basic_analysis(sample, pow_bin=0, autocorr=1, corr_powbin_max=10):
    if autocorr==0:
        ana_func.autocorr_bin(sample, corr_powbin_max),
    mu1 = ana_func.average(sample)
    sigma1 = ana_func.sigma_bootstrap(sample, pow_bin, 300)
    print(f"{mu1} +- {sigma1}")
    ana_func.plot_sampletrend(sample)
    pyp.show()
    return

if __name__=='__main__':
    dati = np.loadtxt("C:\\Users\\e.zippo\\Desktop\\Università\\Tesi\\codice_lattice_su2\\lattice-su2\\su2\\dati_0.55_adj2.4_rnd.dat")
    
    basic_analysis(dati[:,3],0)