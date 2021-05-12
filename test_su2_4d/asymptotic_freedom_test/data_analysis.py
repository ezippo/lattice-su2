import sys
sys.path.append("C:\\Users\\e.zippo\\Desktop\\Università\\Tesi\\codice_lattice_su2\\lattice-su2\\data_analysis")

import numpy as np
from matplotlib import pyplot as pyp
from scipy.optimize import curve_fit

from data_analysis import analysis_functions as ana_func



def basic_analysis(sample, pow_bin=0, autocorr=0, corr_powbin_max=10, bool_plot=0, name=""):
    if autocorr==1:
        ana_func.autocorr_bin(sample, corr_powbin_max),
    mu1 = ana_func.average(sample)
    sigma1 = ana_func.sigma_bootstrap(sample, pow_bin, 100)
    print(f"{name}      {mu1} +- {sigma1}")
    if bool_plot==1:
        ana_func.plot_sampletrend(sample)
        pyp.show()
    return
    
    
def creutz_22(dati, pow_bin=0):
    ''' return creutz ratio W(2,2)W(1,1)/W(1,2)^2 '''
    wloop = dati[:,2]
    wloop_11 = dati[:,5]-1/3
    wloop_10 = dati[:,1]
    
    creutz = ana_func.creutz_ratio(wloop, wloop_11, wloop_10, wloop_10)
    sigma_creutz = ana_func.sigma_bootstrap_creutz(wloop,wloop_11,wloop_10,wloop_10, pow_bin, 500)
    
    return [creutz, sigma_creutz]

def creutz_44(dati, pow_bin=0):
    ''' return creutz ratio W(4,4)W(2,2)/W(4,2)^2 with error '''
    wloop = dati[:,4]
    wloop_11 = dati[:,2]
    wloop_10 = dati[:,3]
    
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
        if n_obs==5:
            mu_sigma_obs[i] = [mu-1/3, sigma]
        else:
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
    beta = np.array([ 6.7, 5., 4., 3.3, 2.8, 2.5,2.35, 2.2, 2., 1.8, 1.66 ])
    beta2 = np.array([ 6.7, 5., 4., 3.3, 2.8, 2.5,2.35, 2.2, 2., 1.8 ])
    
    g0_2 = 4./beta
    g02_2 = 4./beta2
    creutz = np.loadtxt(path_dati+ "\\test_su2_4d\\asymptotic_freedom_test\\wilson_action\\creutz_ratios.txt")
    
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
    

def plot_creutz_ratios_adj(path_file, beta_adj_l, beta_f=0):
    
    g0_2 = 4./(beta_f+8*np.array(beta_adj_l)/3)
    
    creutz = np.loadtxt(path_file+ "creutz_ratios.txt")
    
    #xx0 = np.linspace(0.35,0.8,10)
    #xx1 = np.linspace(0.8,1.2,1000)
    
    pyp.figure()
    #pyp.xlim([0.35,0.9])
    pyp.xlabel(r'$g_0^2$')
    pyp.ylabel(r'1 - $\chi$')
    pyp.title('Creutz ratio')
    pyp.grid()
    #pyp.plot(xx0, 0.04956*xx0, 'k', linewidth=1.2)
    #pyp.plot(xx1, 1-1/xx1, 'k',linewidth=1.2)
    #pyp.plot(xx1, 1-1/xx1**4, 'k',linewidth=1.2)
    pyp.errorbar(g0_2, 1-creutz[:,1], creutz[:,2], fmt='b.', capsize=4, label=r'F($g_0$)' )
    pyp.errorbar(g0_2[:-1], 1-creutz[:-1,3], creutz[:-1,4], fmt='g.', capsize=4, label=r'G($g_0$)')
    pyp.legend(loc='upper left')
    pyp.show()
    

def plot_wloop():
    path_dati = "C:\\Users\\e.zippo\\Desktop\\Università\\Tesi\\codice_lattice_su2\\lattice-su2\\"
    beta = np.array([6.7, 5., 4., 3.3, 2.8, 2.5,2.35, 2.2, 2., 1.8, 1.66 ])
    
    g0_2 = 4./beta
    
    wloop11 = np.loadtxt(path_dati+ "\\test_su2_4d\\asymptotic_freedom_test\\wilson_action\\wloop_1x1.txt")
    wloop21 = np.loadtxt(path_dati+ "\\test_su2_4d\\asymptotic_freedom_test\\wilson_action\\wloop_2x1.txt")
    wloop22 = np.loadtxt(path_dati+ "\\test_su2_4d\\asymptotic_freedom_test\\wilson_action\\wloop_2x2.txt")
    wloop42 = np.loadtxt(path_dati+ "\\test_su2_4d\\asymptotic_freedom_test\\wilson_action\\wloop_4x2.txt")
    wloop44 = np.loadtxt(path_dati+ "\\test_su2_4d\\asymptotic_freedom_test\\wilson_action\\wloop_4x4.txt")
    
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
    
    
def plot_wloop_adj(path_file, beta_adj_l, beta_f=0):
    g0_2 = 4./(beta_f+8*np.array(beta_adj_l)/3)
    
    wloop11 = np.loadtxt(path_file+ "wloop_1x1.txt")
    wloop21 = np.loadtxt(path_file+ "wloop_2x1.txt")
    wloop22 = np.loadtxt(path_file+ "wloop_2x2.txt")
    wloop42 = np.loadtxt(path_file+ "wloop_4x2.txt")
    wloop44 = np.loadtxt(path_file+ "wloop_4x4.txt")
    
    #xx0 = np.linspace(0.43,0.6,1000)
    #xx1 = np.linspace(1.8,3,1000)
    
    pyp.figure()
    pyp.xlabel(r'$g_0^2$')
    pyp.ylabel('W')
    pyp.title('Wilson loop')
    pyp.grid()
    #pyp.plot(xx0, 1-xx0/2, 'k', linewidth=1.2)
    #pyp.plot(xx1, 1/xx1, 'k',linewidth=1.2)
    #pyp.plot(xx1, 1/xx1**2, 'k',linewidth=1.2)
    #pyp.plot(xx1, 1/xx1**4, 'k',linewidth=1.2)
    pyp.errorbar(g0_2, wloop11[:,1], wloop11[:,2], fmt='^',  label='W(1,1)' )
    pyp.errorbar(g0_2, wloop21[:,1], wloop21[:,2], fmt='o', label='W(2,1)' )
    pyp.errorbar(g0_2, wloop22[:,1], wloop22[:,2], fmt='p', label='W(2,2)' )
    pyp.errorbar(g0_2, wloop42[:,1], wloop42[:,2], fmt='s',  label='W(4,2)' )
    pyp.errorbar(g0_2, wloop44[:,1], wloop44[:,2], fmt='v',  label='W(4,4)' )
    pyp.legend()
    pyp.show()
    

def plot_asyfree(path_file, beta_l):
    
    g0_2 = 4./np.array(beta_l)
    gr_2 = 1./ (1./g0_2 - 11*np.log(2)/(12*np.pi*np.pi) )
    
    creutz = np.loadtxt(path_file + "creutz_ratios.txt")
    xx0 = np.linspace(0,1,10)
    
    pyp.figure()
    pyp.xlabel(r'$g_0^2$')
    pyp.ylabel(r'1 - $\chi$')
    pyp.title('Asymptotic freedom')
    pyp.grid()
    pyp.plot(xx0, 0.04956*xx0, 'k', linewidth=1.2)
    pyp.errorbar(g0_2, 1-creutz[:,1], creutz[:,2], fmt='b.', capsize=4, label=r'F($g_0$)' )
    pyp.errorbar(gr_2[:-1], 1-creutz[:-1,3], creutz[:-1,4], fmt='g.', capsize=4, label=r'G($(\frac{1}{g_0^2} - \frac{11 ln2}{12 \pi^2})^{-1/2}$)')
    pyp.legend(loc='upper left')
    pyp.show()
    
    

def plot_asyfree_adj(path_file, beta_adj_l, beta_f=0):
    
    g0_2 = 4./(beta_f+8*np.array(beta_adj_l)/3)
    gr_2 = 1./ (1./g0_2 - 11*np.log(2)/(12*np.pi*np.pi) )
    
    creutz = np.loadtxt(path_file+ "creutz_ratios.txt")
    #xx0 = np.linspace(0,1,10)
    
    pyp.figure()
    #pyp.xlim([0.2,0.9])
    pyp.xlabel(r'$g_0^2$')
    pyp.ylabel(r'1 - $\chi$')
    pyp.title('Asymptotic freedom')
    pyp.grid()
    #pyp.plot(xx0, 0.04956*xx0, 'k', linewidth=1.2)
    pyp.errorbar(g0_2, 1-creutz[:,1], creutz[:,2], fmt='b.', capsize=4, label=r'F($g_0$)' )
    pyp.errorbar(gr_2[:-1], 1-creutz[:-1,3], creutz[:-1,4], fmt='g.', capsize=4, label=r'G($(\frac{1}{g_0^2} - \frac{11 ln2}{12 \pi^2})^{-1/2}$)')
    pyp.legend(loc='upper left')
    pyp.show()
    
    
def plot_gr_vs_g0():
    path_dati = "C:\\Users\\e.zippo\\Desktop\\Università\\Tesi\\codice_lattice_su2\\lattice-su2\\"
    beta = np.array([ 6.7, 5., 4., 3.3, 2.8, 2.5,2.35, 2.2, 2., 1.8, 1.66 ])
    
    invg0_2 = beta/4.
    
    creutz = np.loadtxt(path_dati+ "\\test_su2_4d\\asymptotic_freedom_test\\wilson_action\\creutz_ratios.txt")
    invgr_2 = 0.04956/(1-creutz[:,1])
    d_invgr2 = 0.04956*creutz[:,2]/(1-creutz[:,1])**2
    print(invg0_2)
    
    def retta(x,a,b):
        return a*x + b
    
    n=4
    popt, pcov = curve_fit(retta, invg0_2[:n], invgr_2[:n], sigma=d_invgr2[:n])
    print(popt)
    print(np.sqrt(pcov.diagonal()))
    chi2 = np.sqrt( ( ((retta(invg0_2[:n], *popt)-invgr_2[:n])/d_invgr2[:n])**2 ).sum() )
    print(chi2)
    xx0 = np.linspace(0,0.6,1000)
    xx1 = np.linspace(0.3,2,10)
    
    pyp.figure()
    pyp.xlabel(r'$1/g_0^2$')
    pyp.ylabel(r'1/$g^2(2a)$')
 #   pyp.title('Asymptotic freedom')
    pyp.grid()
    pyp.plot(xx1, retta(xx1, popt[0], popt[1]), 'k', linewidth=1.1)
    pyp.plot(xx0, 0.04956/(1-xx0), 'k', linewidth=1.1)
    pyp.errorbar(invg0_2, invgr_2, d_invgr2, fmt='r.', capsize=4 )
 #   pyp.errorbar(gr_2[:-1], 1-creutz[:-1,3], creutz[:-1,4], fmt='g.', capsize=4, label=r'G($(\frac{1}{g_0^2} - \frac{11 ln2}{12 \pi^2})^{-1/2}$)')
  #  pyp.legend(loc='upper left')
    pyp.show()
    

def plot_gr_vs_g0_adj():
    path_dati = "C:\\Users\\e.zippo\\Desktop\\Università\\Tesi\\codice_lattice_su2\\lattice-su2\\"
    beta_adj = np.array([3., 2.2, 1.8, 1.7, 1.65, 1.6, 1.57, 1.55, 1.53, 1.5, 1.4 ])
    beta_fund = 1.0
    
    invg0_2 = (beta_fund + 8*beta_adj/3)/4.
    
    creutz = np.loadtxt(path_dati+ "\\test_su2_4d\\asymptotic_freedom_test\\fund_plus_adj_action\\creutz_ratios.txt")
    invgr_2 = 0.04956/(1-creutz[:,1])
    d_invgr2 = 0.04956*creutz[:,2]/(1-creutz[:,1])**2
    print(invg0_2)
    
    def retta(x,a,b):
        return a*x + b
    
    n=4
    popt, pcov = curve_fit(retta, invg0_2[:n], invgr_2[:n], sigma=d_invgr2[:n])
    print(popt)
    print(np.sqrt(pcov.diagonal()))
    chi2 = np.sqrt( ( ((retta(invg0_2[:n], *popt)-invgr_2[:n])/d_invgr2[:n])**2 ).sum() )
    print(chi2)
    #xx0 = np.linspace(1.1,1.4,1000)
    xx1 = np.linspace(0.8,2.3,10)
    
    pyp.figure()
    pyp.xlabel(r'$1/g_0^2$')
    pyp.ylabel(r'1/$g^2(2a)$')
 #   pyp.title('Asymptotic freedom')
    pyp.grid()
    pyp.plot(xx1, retta(xx1, popt[0], popt[1]), 'k', linewidth=1.1)
    #pyp.plot(xx0, 0.04956/(1-xx0), 'k', linewidth=1.1)
    pyp.errorbar(invg0_2, invgr_2, d_invgr2, fmt='r.', capsize=4 )
 #   pyp.errorbar(gr_2[:-1], 1-creutz[:-1,3], creutz[:-1,4], fmt='g.', capsize=4, label=r'G($(\frac{1}{g_0^2} - \frac{11 ln2}{12 \pi^2})^{-1/2}$)')
  #  pyp.legend(loc='upper left')
    pyp.show()
    

if __name__=='__main__':
    path_dati = "C:\\Users\\e.zippo\\Desktop\\Università\\Tesi\\codice_lattice_su2\\lattice-su2\\"
    #dati = np.loadtxt(path_dati + "su2\\dati_adj2.48.dat")
    #basic_analysis(dati[:,2],0,1,10,1)
    beta = [3., 2.7, 2.6, 2.55, 2.5, 2.48, 2.465, 2.45, 2.4]
    powbin = [5., 6., 6., 7., 6., 6.,  6., 6., 5.]
    file_l = ['3' , '2.7', '2.6', '2.55','2.5', '2.48', '2.465', '2.45', '2.4']
    
    #creutz_savetxt("creutz_ratios.txt", path_dati+"su2\\dati_adj", file_l, beta, powbin)
    #plot_creutz_ratios_adj(path_dati+"test_su2_4d\\asymptotic_freedom_test\\adj_action\\", beta)
    #plot_wloop_adj(path_dati+"test_su2_4d\\asymptotic_freedom_test\\adj_action\\", beta)
    #creutz_savetxt("creutz_ratios.txt",path_dati+"su2\\dati_adj", file_l, beta, powbin) 
    #plot_asyfree_adj(path_dati+"test_su2_4d\\asymptotic_freedom_test\\adj_action\\", beta)
    '''
    path_dati = "C:\\Users\\e.zippo\\Desktop\\Università\\Tesi\\codice_lattice_su2\\lattice-su2\\"
    dati = np.loadtxt(path_dati + "su2\\dati_1_adj1.51.dat")
    beta_adj = [3., 2.2, 1.8, 1.7, 1.65, 1.6, 1.57, 1.55, 1.53, 1.5, 1.4]
    file_list = ['3', '2.2', '1.8', '1.7', '1.65', '1.6', '1.57', '1.55', '1.53', '1.5', '1.4']
    powbin = np.loadtxt(path_dati+"test_su2_4d\\asymptotic_freedom_test\\fund_plus_adj_action\\pow_bin_list.txt")
    
    plot_gr_vs_g0_adj()
    '''