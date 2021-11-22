import subprocess
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

palette = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

def main():
    min_N = input("Enter the minimum number of rows: ")
    max_N = input("Enter the maximum number of rows: ")
    # * Matrix product executable file is assumed to be present
    proc = subprocess.Popen("./matrix_product",stdin = subprocess.PIPE,stdout = subprocess.PIPE)
    proc.stdin.write(f'100,{min_N},{max_N}'.encode('UTF-8'))
    proc.stdin.close()
    # TODO: adapt for unsuccessful timing scenario

    '''Importing data from performance.dat'''
    performance = np.genfromtxt('performance.dat')
    fig, ax = plt.subplots()
    x = performance[:,0]
    for i in range(1,4):
        '''Leading order term fitting'''
        f_theo = lambda x, a, b, c: a*x**b + c
        y = performance[:,i]
        (a,b,c), _ = curve_fit(f_theo, x, y, maxfev = 2000)
        f = lambda x: f_theo(x,a,b,c)
        _, _ax = plt.subplots()
        _ax.scatter(x,y, color = palette[i])
        _ax.plot(x,f(x))
        _ax.set_title(f'Empirical complexity method #{i}')
        _ax.set_xlabel('input size N')
        _ax.set_ylabel('computation time (seconds)')
        power = f'{round(b,2)}'
        _ax.text(0, max(y),r'Complexity: $\mathcal{O}\left(n^{'+f'{power}'+r'}\right)$', fontsize=10)
        _.savefig(f'performance_method{i}.jpg')
        ax.scatter(x,y, color = palette[i])
    ax.legend([f'method {i}' for i in range(1,4)])
    ax.set_xlabel('input size N')
    ax.set_ylabel('computation time (seconds)')
    ax.set_title('Performance comparison')
    fig.savefig(f'performances.jpg')
    print('Successfully generated plots')

if __name__ == '__main__':
    main()
