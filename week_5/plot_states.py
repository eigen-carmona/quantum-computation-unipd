import subprocess
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import time


def main():
    '''Executing fortran program'''
    L = 10
    dts = [0.1,0.05,0.01]
    dx = 2*L/(2**12)
    Ts = [4,8,16]
    states = [0,1,2,3,4]
    energy = lambda n: (1/2+n)
    x = np.arange(-L+dx,L+dx,dx)
    for dt in dts:
        for T in Ts:
            steps = int(T/dt)
            t = np.arange(0,(steps+1)*dt,dt)
            psi_ev = np.empty((len(x),steps+1,len(states)), dtype = complex)
            for n in states:
                proc = subprocess.Popen(
                    "./time_dependent",
                    stdin = subprocess.PIPE,
                    stdout = subprocess.PIPE
                    )
                out = proc.communicate(f'{n},{steps},{T}'.encode('UTF-8'))[0]
                # Load the resulting data
                psi_ev[:,:,n] = np.loadtxt(
                    f'data/state_{n}_evo_T{T}_{steps}.dat',
                    dtype=complex,
                    delimiter=',',
                    converters = {-1: lambda s: complex(0)}
                    ).T
            fig, ax = plt.subplots()
            lines = [None]*len(states)
            line, = ax.plot(x,(x-t[0]/T)**2/2)
            pdf = np.abs(psi_ev)**2
            _fig, _ax = plt.subplots()
            for n in states:
                avg_x = dx*np.array([x]).T*pdf[:,:,n]
                avg_x = avg_x.sum(axis = 0)
                _ax.plot(t,avg_x)
                lines[n], = ax.plot(x, pdf[:,0,n])
            _ax.set_title(r'Average position for $T='+f'{T}'+r'$ and $dt='+f'{dt}'+r'$')
            _ax.set_xlabel(r'$t$')
            _ax.set_ylabel(r'$\langle x\rangle$')
            _fig.savefig(f'plots/states_T{T}_{steps}.jpg')
            ymax = 1.5*max(pdf[:,0,4]+energy(4))
            ax.set_ylim(ymin = 0, ymax = ymax)
            ax.legend(['potential']+[r'$n ='+f'{n}'+r'$' for n in states])
            ax.set_xlabel('x')
            ax.set_ylabel(r'$|\psi|^2$')
            ax.set_title(r'Time evolution for $T='+f'{T}'+r'$ and $dt='+f'{dt}'+r'$')
            def animate(i):
                line.set_ydata((x-t[i]/T)**2/2)
                for n in states:
                    lines[n].set_ydata(pdf[:,i,n]+energy(n))
                return tuple([line]+lines)
            interval =int(2000/steps)
            ani = animation.FuncAnimation(fig, animate, interval = interval, frames = steps)
            ani.save(f'animations/states_T{T}_{steps}.gif')
            plt.close()

if __name__ == '__main__':
    main()
