import subprocess
import matplotlib.pyplot as plt
import numpy as np
import time


def main():
    Js = np.arange(0,3.2,0.2)
    spins_n = [i for i in range(3,12)]
    m_sizes = [4,8,12,16]
    k_levels = [min(2**n,5) for n in spins_n]
    timing = []
    '''Executing fortran program'''
    # Plot separately the first 4 energy levels as a function of J for the 8 spin model
    ## Load the exact 8 particle Ising model energies
    energies = np.empty((len(Js),2**8))
    _energies = np.empty((len(Js),max(m_sizes),len(m_sizes)))
    for i,j in enumerate(Js):
        energies[i,:] = np.loadtxt(f'../week_7/data/energies_8_spins_lambda_{"%.2f"%j}.dat').T
        ## Generate the datasets for the different values of J and the different m
        for k,m in enumerate(m_sizes):
            #proc = subprocess.Popen(
            #    "./real_group_ising.exe",
            #    stdin = subprocess.PIPE,
            #    stdout = subprocess.PIPE
            #    )
            ## (initial system size, states in approximation, lambda value, iterations)
            #out = proc.communicate(f'4,{m},{j},1'.encode('UTF-8'))[0]
            ### Load the rg estimated energies
            _energies[i,:m,k] = np.loadtxt(f'data/{m}_energies_8_spins_lambda_{"%.2f"%j}.dat').T
    ## Plot for each energy level
    for n_state in range(4):
        legend = []
        _fig, _ax = plt.subplots()
        for k,m in enumerate(m_sizes):
            ## Measure error
            rel_error = np.abs((energies[:,n_state]-_energies[:,n_state,k])/energies[:,n_state])
            _ax.plot(Js[1:],np.log(rel_error[1:]))
            rel_error = np.average(rel_error)
            legend.append(r'$m = '+f'{m}'+r'$, $\log{(<err_{rel}>)}: '+f'{"%.2f"%np.log(rel_error)}'+r'$')
        _ax.legend(legend)
        _ax.set_title(r'Energy estimations for $n\_state='+f'{n_state+1}'+r'$')
        _ax.set_ylabel(r'$\log($relative error$)$')
        _ax.set_xlabel(r'$\lambda$')
        _fig.savefig(f'plots/8_energy_error_{n_state+1}.jpg')
    # Plot the convergence
    fig, ax = plt.subplots()
    # Plot the energy convergence
    _energies = np.empty((len(Js),max(m_sizes),len(m_sizes)))
    _fig, _ax = plt.subplots()
    legend = []
    for k,m in enumerate(m_sizes):
        iterations = []
        for i,j in enumerate(Js):
            proc = subprocess.Popen(
                "./real_group_ising.exe",
                stdin = subprocess.PIPE,
                stdout = subprocess.PIPE
                )
            # (initial system size, states in approximation, lambda value, iterations)
            out = proc.communicate(f'8,{m},{j},0'.encode('UTF-8'))[0]
            # Program response is the required iterations for convergence
            total_iters = int(out.decode('UTF-8').split('\n')[-2])
            iterations.append(total_iters)
            n_spins = 8*2**total_iters
            _energies[i,:m,k] = np.loadtxt(f'data/{m}_energies_{n_spins}_spins_lambda_{"%.2f"%j}.dat').T/n_spins
        ax.plot(Js,iterations)
        _ax.plot(Js,_energies[:,0,k])
        legend.append(r'$m='+f'{m}'+r'$')
    ax.legend(legend)
    ax.set_title(r'Convergence for $\epsilon=10^{-10}$')
    ax.set_xlabel(r'$\lambda$')
    ax.set_ylabel(r'Minimum iterations')
    fig.savefig('plots/convergence_iterations.jpg')

    _ax.legend(legend)
    _ax.set_title(r'Energy density convergence')
    _ax.set_xlabel(r'$\lambda$')
    _ax.set_ylabel(r'$e\ (E/N)$')
    _fig.savefig('plots/energy_convergence.jpg')


if __name__ == '__main__':
    main()
