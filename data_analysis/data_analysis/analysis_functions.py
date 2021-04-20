import time

import numpy as np
from matplotlib import pyplot as pyp
from scipy.optimize import curve_fit


def average(sample, n_term=0):
    ''' sample average 
        parameter n_term: first n_term data are discarded
        '''
    return (sample[n_term:]).sum()/(len(sample)-n_term)
    
def sigma_simple(sample, n_term=0):
    ''' sample standard deviation with zero covariance
        parameter n_term: first n_term data are discarded
        '''
    mu = average(sample, n_term)
    N = len(sample) - n_term
    return np.sqrt(((sample[n_term:] - mu)**2).sum()/(N*(N-1)))
    
def sigma_blocking(sample, pow_bin=0,  n_term=0):
    ''' sample standard deviation with covariance estimate with data blocking technique 
        parameter n_term: first n_term data are discarded
        parameter pow_bin: block size = 2^pow_bin
        '''
    dim_bin = 2**pow_bin
    n_bin = int(( len(sample[n_term:]) )/ dim_bin)    
    N = n_bin * dim_bin
    
    mu = average(sample[n_term:n_term+N])
    sample_split = np.array_split(sample[n_term:n_term+N],n_bin)
    sum_block = (np.array([ ( average(sample[n_term+i*dim_bin : n_term+(i+1)*dim_bin]) - mu )**2 for i in range(n_bin) ]) ).sum()
    
    return np.sqrt(sum_block)*dim_bin/N
    

def autocorr_bin(sample, pow_bin_max=10, n_term=0):
    ''' plot of sigma_blocking for different pow_bin:
            correlation lenght estimate
            '''
    pow_bin = [x for x in range(pow_bin_max)]
    sigma = np.array([sigma_blocking(sample, x, n_term) for x in pow_bin])
    sigma = sigma/sigma[0]
    pyp.figure(10)
    #pyp.yscale('log')
    pyp.grid(True)
    pyp.plot(pow_bin, sigma)
    pyp.show()


def sigma_bootstrap(sample, pow_bin=0, n_resample=100, function='average', n_term=0):
    ''' estimator standard deviation with covariance estimate with bootstrap technique
        parameter pow_bin: block size = 2^pow_bin
        parameter n_resample: number of resamples from sample
        parameter function: function of the estimator
        parameter n_term: first n_term data are discarded
        '''
    dim_bin = 2**pow_bin
    n_bin = int(( len(sample[n_term:]) )/ dim_bin)    
    N = n_bin * dim_bin
    
    sample_split = np.reshape( sample[n_term: n_term+N], (n_bin, dim_bin) )
    rnd_matrix = np.random.randint(n_bin, size=(n_resample,n_bin))
    resample = np.array([ np.reshape( np.array([ sample_split[i] for i in rnd_matrix[k] ]), N ) for k in range(n_resample) ])
    
    if function=='average':
        mu = np.array([ average(resample[i]) for i in range(n_resample) ])
    else:
        print("ERROR: wrong function in sigma_bootstrap()")
    
    return np.sqrt(((mu - average(mu))**2).sum()/n_resample)


def plot_sampletrend(sample, x_label="x",y_label="y",title="", figure=100 ):
    ''' plot the measures evolution at each MonteCarlo step '''
    n = [i for i in range (len(sample))]
    pyp.figure(figure)
    pyp.plot(n,sample)
    pyp.grid(True)
    pyp.xlabel(x_label)
    pyp.ylabel(y_label)
    pyp.title(title)
    

def plot_errorbar(x,y,dy, x_label="x", y_label="y", title="", dx=0, fmt_=".", figure=200 ):
    ''' plot with error '''
    pyp.figure(figure)
    pyp.errorbar(x,y,dy,dx,fmt=fmt_)
    pyp.xlabel(x_label)
    pyp.ylabel(y_label)
    pyp.grid(True)
    pyp.title(title)


#if __name__=='__main__':
    