import subprocess
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import time


def main():
    '''Starts the fortran matrix_product and enters the user specified values.
    Generates plots for the performance measurements
    and fits them to estimate algorithm complexity.'''
    samples = input("How many matrices should be generated? ")
    # ! Matrix product executable file is assumed to be present
    '''Executing fortran program'''
    for i in range(int(samples)):
        proc = subprocess.Popen("../eigenproblem/./eigenproblem",stdin = subprocess.PIPE,stdout = subprocess.PIPE)
        out = proc.communicate(f'1000,"sample{i}.dat"'.encode('UTF-8'))[0]
    print(f'generated {samples} datasets...')

    '''Retrieving normalized spacings'''
    v = np.genfromtxt('sample0.dat')
    for i in range(1,int(samples)):
        # Horizontally stacking spacing samples
        v = np.hstack((v,np.genfromtxt(f'sample{i}.dat')))
    # numpy histogram to fit
    y, x = np.histogram(v, bins = 500)
    p_s_th = lambda s, A, a, B, b: A*(s**a)*np.exp(-B*s**b)
    (A,a,B,b), _ = curve_fit(p_s_th, x[:-1], y)

    # plotting histogram and fitted curve
    plt.hist(v,bins = 500)
    p_s = lambda s: p_s_th(s,A,a,B,b)
    plt.plot(x[:-1],p_s(x[:-1]))
    plt.title('Spacings frequencies')
    plt.xlabel('s')
    plt.ylabel('N')
    plt.text(max(x[:-int(len(x)/2)]), max(y),\
        r'$'+f'{round(A,2)}'+\
        r's^{'+f'{round(a,2)}'+\
        r'}\exp\left(-'+ f'{round(B,2)}'+\
        r's^{'+f'{round(b,2)}' +r'}\right)$',\
        fontsize=10)
    plt.savefig('histogram.jpg')
    

if __name__ == '__main__':
    main()
