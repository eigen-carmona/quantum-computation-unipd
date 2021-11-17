import subprocess, numpy, matplotlib

def main():
    min_N = input("Enter the minimum number of rows: ")
    max_N = input("Enter the maximum number of rows: ")
    # * Matrix product executable file is assumed to be present
    proc = subprocess.Popen("./matrix_product",stdin = subprocess.PIPE,stdout = subprocess.PIPE)
    proc.stdin.write(f'100,{min_N},{max_N}'.encode('UTF-8'))
    proc.stdin.close()

if __name__ == '__main__':
    main()
