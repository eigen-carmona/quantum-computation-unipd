import subprocess
import matplotlib.pyplot as plt
import numpy as np

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
    for i in range(1,4):
        plot = plt.plot(performance[:,0],performance[:,i])
        plt.savefig(f'performance_method{i}.jpg')
    print('Successfully generated plots')

if __name__ == '__main__':
    main()
