import sys
sys.path.append("C:\\Users\\e.zippo\\Desktop\\Università\\Tesi\\codice_lattice_su2\\lattice-su2\\data_analysis")

import numpy as np
from matplotlib import pyplot as pyp
from scipy.optimize import curve_fit

from data_analysis import analysis_functions as ana_func



def basic_analysis(sample, pow_bin=0, autocorr=1, corr_powbin_max=10, name=""):
    if autocorr==0:
        ana_func.autocorr_bin(sample, corr_powbin_max),
    mu1 = ana_func.average(sample)
    sigma1 = ana_func.sigma_bootstrap(sample, pow_bin, 300)
    print(f"{name}      {mu1} +- {sigma1}")
    ana_func.plot_sampletrend(sample)
    pyp.show()
    return
    
    
def creutz_22(dati, pow_bin=0):
    ''' return creutz ratio W(2,2)W(1,1)/W(1,2)^2 '''
    wloop = dati[:,1]
    wloop_11 = dati[:,2]
    wloop_10 = dati[:,3]
    
    creutz = ana_func.creutz_ratio(wloop, wloop_11, wloop_10, wloop_10)
    sigma_creutz = ana_func.sigma_bootstrap_creutz(wloop,wloop_11,wloop_10,wloop_10, pow_bin, 500)
    
    return [creutz, sigma_creutz]

def creutz_44(dati, pow_bin=0):
    ''' return creutz ratio W(4,4)W(2,2)/W(4,2)^2 with error '''
    wloop = dati[:,4]
    wloop_11 = dati[:,1]
    wloop_10 = dati[:,5]
    
    creutz = ana_func.creutz_ratio(wloop, wloop_11, wloop_10, wloop_10)
    sigma_creutz = ana_func.sigma_bootstrap_creutz(wloop,wloop_11,wloop_10,wloop_10, pow_bin, 500)
    
    return [creutz, sigma_creutz]
    
    
def mu_and_sigma_data(path_file, file_name, pow_bin_l, n_obs):
    ''' returns an array of [average,sigma] of the observable dati[:,n_obs] for each data file in list file_name '''
    if type(path_file) is not type(''):
        print("ERROR: parameter path_file in mu_and_sigma_dati() must be a string")
        return 0
    n_dati = len(file_name)
    mu_sigma_obs = np.zeros((n_dati,2))
    for i in range(n_dati):
        dati = np.loadtxt(path_file + file_name[i] +".dat")
        mu = ana_func.average(dati[:,n_obs])
        sigma= ana_func.sigma_bootstrap(dati[:,n_obs], pow_bin_l[i], 500)
        mu_sigma_obs[i] = [mu, sigma]
    
    return mu_sigma_obs
    
    
def wloop_savetxt(txt_name, path_file, data_name_list, beta_list, pow_bin_l, n_obs, fmt='%.6f'):
    ''' save on a txt an array of [beta,average,sigma] of the observable dati[:,n_obs] for each beta in beta_list '''
    if type(path_file) is not type('') or type(txt_name) is not type(''):
        print("ERROR: parameter path_file in mu_and_sigma_dati() must be a string")
        return 0
    wloop = mu_and_sigma_data(path_file, data_name_list, pow_bin_l, n_obs)
    array_to_record = np.array( [ [beta_list[i], wloop[i][0], wloop[i][1]] for i in range(len(beta_list)) ] )
    np.savetxt(txt_name, array_to_record, fmt)
    
    
def creutz_from_data(path_file, file_name, pow_bin_l, n_cr=0):
    ''' returns an array of [creutz,sigma] for each data file in list file_name
        creutz_22 if n_cr=0, creutz_44 if n_cr=1'''
    if type(path_file) is not type(''):
        print("ERROR: parameter path_file in mu_and_sigma_dati() must be a string")
        return 0
    n_dati = len(file_name)
    creutz_arr = np.zeros((n_dati,2))
    for i in range(n_dati):
        dati = np.loadtxt(path_file + file_name[i] +".dat")
        if n_cr==0:
            creutz_arr[i] = creutz_22(dati,pow_bin_l[i])
        elif n_cr==1:
            creutz_arr[i] = creutz_44(dati,pow_bin_l[i])
        else:
            print('ERROR: parameter n_cr in creutz_from_data() must be 0 or 1 ')
    
    return creutz_arr
    

def creutz_savetxt(txt_name, path_file, data_name_list, beta_list, pow_bin_l, fmt='%.6f'):
    ''' save on a txt an array of [beta,creutz22,sigma22,creutz44,sigma44] for each beta in beta_list '''
    if type(path_file) is not type('') or type(txt_name) is not type(''):
        print("ERROR: parameter path_file in mu_and_sigma_dati() must be a string")
        return 0
    creutz = creutz_from_data(path_file, data_name_list, pow_bin_l, 0)
    creutz2 = creutz_from_data(path_file, data_name_list, pow_bin_l, 1)
    array_to_record = np.array([ [beta_list[i], creutz[i][0], creutz[i][1], creutz2[i][0], creutz2[i][1]] for i in range(len(beta_list)) ])
    np.savetxt(txt_name, array_to_record, fmt)
    

def plot_creutz_ratios():
    path_dati = "C:\\Users\\e.zippo\\Desktop\\Università\\Tesi\\codice_lattice_su2\\lattice-su2\\"
    file_name_l = ['6.7', '5', '4', '3.3', '2.8', '2.5', '2.35', '2.2', '2', '1.8', '1.66']
    beta = np.array([ 6.7, 5., 4., 3.3, 2.8, 2.5,2.35, 2.2, 2., 1.8, 1.66 ])
    beta2 = np.array([ 6.7, 5., 4., 3.3, 2.8, 2.5,2.35, 2.2, 2., 1.8 ])
    powbin = np.loadtxt(path_dati+"test_su2_4d\\asymptotic_freedom_test\\pow_bin_list.txt")
    
    g0_2 = 4./beta
    g02_2 = 4./beta2
    creutz = np.loadtxt(path_dati+ "\\test_su2_4d\\asymptotic_freedom_test\\creutz_ratios.txt")
    
    xx0 = np.linspace(0,1.5,10)
    xx1 = np.linspace(2,3,1000)
    
    pyp.figure()
    pyp.xlabel(r'$g_0^2$')
    pyp.ylabel(r'1 - $\chi$')
    pyp.title('Creutz ratio')
    pyp.grid()
    pyp.plot(xx0, 0.04956*xx0, 'k', linewidth=1.2)
    pyp.plot(xx1, 1-1/xx1, 'k',linewidth=1.2)
    pyp.plot(xx1, 1-1/xx1**4, 'k',linewidth=1.2)
    pyp.errorbar(g0_2, 1-creutz[:,1], creutz[:,2], fmt='b.', capsize=4, label=r'F($g_0$)' )
    pyp.errorbar(g02_2, 1-creutz[:-1,3], creutz[:-1,4], fmt='g^', capsize=4, label=r'G($g_0$)')
    pyp.legend(loc='upper left')
    pyp.show()
    

def plot_wloop():
    path_dati = "C:\\Users\\e.zippo\\Desktop\\Università\\Tesi\\codice_lattice_su2\\lattice-su2\\"
    beta = np.array([6.7, 5., 4., 3.3, 2.8, 2.5,2.35, 2.2, 2., 1.8, 1.66 ])
    
    g0_2 = 4./beta
    
    wloop11 = np.loadtxt(path_dati+ "\\test_su2_4d\\asymptotic_freedom_test\\wloop_1x1.txt")
    wloop21 = np.loadtxt(path_dati+ "\\test_su2_4d\\asymptotic_freedom_test\\wloop_2x1.txt")
    wloop22 = np.loadtxt(path_dati+ "\\test_su2_4d\\asymptotic_freedom_test\\wloop_2x2.txt")
    wloop42 = np.loadtxt(path_dati+ "\\test_su2_4d\\asymptotic_freedom_test\\wloop_4x2.txt")
    wloop44 = np.loadtxt(path_dati+ "\\test_su2_4d\\asymptotic_freedom_test\\wloop_4x4.txt")
    
    xx0 = np.linspace(0,1,1000)
    xx1 = np.linspace(1.8,3,1000)
    
    pyp.figure()
    pyp.xlabel(r'$g_0^2$')
    pyp.ylabel('W')
    pyp.title('Wilson loop')
    pyp.grid()
    pyp.plot(xx0, 1-3*xx0/16, 'k', linewidth=1.2)
    pyp.plot(xx1, 1/xx1, 'k',linewidth=1.2)
    pyp.plot(xx1, 1/xx1**2, 'k',linewidth=1.2)
    pyp.plot(xx1, 1/xx1**4, 'k',linewidth=1.2)
    pyp.errorbar(g0_2, wloop11[:,1], wloop11[:,2], fmt='^',  label='W(1,1)' )
    pyp.errorbar(g0_2, wloop21[:,1], wloop21[:,2], fmt='o', label='W(2,1)' )
    pyp.errorbar(g0_2, wloop22[:,1], wloop22[:,2], fmt='p', label='W(2,2)' )
    pyp.errorbar(g0_2, wloop42[:,1], wloop42[:,2], fmt='s',  label='W(4,2)' )
    pyp.errorbar(g0_2, wloop44[:,1], wloop44[:,2], fmt='v',  label='W(4,4)' )
    pyp.legend()
    pyp.show()
    

if __name__=='__main__':
    plot_wloop()
    '''
    path_dati = "C:\\Users\\e.zippo\\Desktop\\Università\\Tesi\\codice_lattice_su2\\lattice-su2\\"
    file_list = ['6.7', '5', '4', '3.3', '2.8', '2.5', '2.35', '2.2', '2', '1.8', '1.66']
    beta = np.array([6.7, 5., 4., 3.3, 2.8, 2.5, 2.35, 2.2, 2, 1.8, 1.66])
    g_02 = 4./beta
    
    powbin = [4, 6, 6,	6,	7, 9, 8, 7,	5, 4, 3]

    creutz_savetxt("creutz_ratios.txt",path_dati+"\\su2\\dati_", file_list, beta, powbin)
    '''