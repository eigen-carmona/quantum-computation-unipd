import subprocess
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import time


def main():
    '''Executing fortran program'''
    sizes = [10000,20000,30000,40000,50000]
    ho = lambda x: x**2/2
    heaviside = lambda x: 1 if x > 0 else 0
    qs = lambda x: np.array([heaviside(i) for i in x])
    qb = lambda x: np.array([-heaviside(i-5)+heaviside(i+5) for i in x])
    qw = lambda x: -qb(x)
    bx = lambda x: 0*qb(x)
    potentials = {'ho': ho, 'qs':qs, 'qb':qb, 'qw': qw, 'bx':bx}
    states = range(100)
    ho_energies = [n+0.5 for n in states]
    x_0 = -10
    x_n = 10
    for pot,potential in potentials.items():
        fig, ax = plt.subplots()
        _fig, _ax = plt.subplots()
        __fig, __ax = plt.subplots()
        for size in sizes:
            proc = subprocess.Popen("./time_independent",stdin = subprocess.PIPE,stdout = subprocess.PIPE)
            out = proc.communicate(f'"{pot}",{size},100'.encode('UTF-8'))[0]
            energies = np.genfromtxt(f'energies_{pot}_{size}.dat')
            psi_states = np.genfromtxt(f'psi_states_{pot}_{size}.dat').T
            dx = (x_n-x_0)/size
            x = np.linspace(x_0+dx,x_n-dx,size-1)
            _ax.plot(x,energies[4]+np.abs(psi_states[:,4])**2/sum(dx*np.abs(psi_states[:,4])**2))
            if pot == 'ho':
                __ax.plot(range(100),energies)
        if pot == 'ho':
            __ax.plot(range(100),ho_energies)
            __ax.set_title(f'{pot} potential, energy predictions')
            __ax.set_xlabel('state n')
            __ax.set_ylabel('Energy')
            __ax.legend([f'N = {size}' for size in sizes]+['Analytic expression'])
            __fig.savefig('Energy fit ho.jpg')
        # plot the third state in different grid sizes
        _ax.plot(x,potential(x))
        _ax.set_title(f'{pot} potential, state n=3')
        _ax.set_xlabel('x')
        _ax.set_ylabel(r'$|\psi|^2$')
        _ax.legend([f'N ={i}' for i in sizes]+['V(x)'])
        _fig.savefig(f'psi3_{pot}.jpg')    
        # Plots: plot the first five states for each potential at grid size 25000
        for i in range(5):
            ax.plot(x,energies[i]+np.abs(psi_states[:,i])**2/sum(dx*np.abs(psi_states[:,i])**2))
        ax.plot(x,potential(x))
        ax.set_title(f'{pot} potential')
        ax.set_xlabel('x')
        ax.set_ylabel(r'$|\psi|^2$')
        ax.legend([f'state {i}' for i in range(5)]+['V(x)'])
        fig.savefig(f'wavefunction_{pot}.jpg')

if __name__ == '__main__':
    main()
