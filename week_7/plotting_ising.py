import subprocess
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import time


def main():
    Js = np.arange(0,3.2,0.2)
    spins_n = [i for i in range(3,12)]
    k_levels = [min(2**n,5) for n in spins_n]
    '''Executing fortran program'''
    # Plot the first k energy levels as a function of J
    for n,k in zip(spins_n,k_levels):
        # Empty array with energy levels for different Js
        energies = np.empty((len(Js),2**n))
        for i,j in enumerate(Js):
            proc = subprocess.Popen(
                "./ising_model.exe",
                stdin = subprocess.PIPE,
                stdout = subprocess.PIPE
                )
            out = proc.communicate(f'{n},{j}'.encode('UTF-8'))[0]
            # Load the resulting data
            energies[i,:] = np.loadtxt(f'data/energies_{n}_spins_lambda_{"%.2f"%j}.dat').T

        fig, ax = plt.subplots()
        for i in range(k):
            ax.plot(Js,energies[:,i])
        ax.set_title(r'Energy levels for $'+f'{n}'+r'$ spins')
        ax.set_ylabel(r'$E$')
        ax.set_xlabel(r'$\lambda$')
        ax.legend([r'$E_'+f'{i}'+r'$' for i in range(k)])
        fig.savefig(f'plots/energy_{n}_spins.jpg')
        plt.close()

if __name__ == '__main__':
    main()
