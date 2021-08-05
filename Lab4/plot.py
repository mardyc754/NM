import numpy as np
import matplotlib.pyplot as plt

def plot(filename):
    fig = plt.figure()

    data = np.loadtxt(filename)

    alfa = data[:, 0]	
    vals = data[:, 1:]
    
    if filename == "out1.txt":
        for i in np.arange(6):
            plt.plot(alfa, vals[:,i], marker='.', label=r"$\omega_{{ {0} }} = \sqrt{{ \lambda_{{ {0} }} }}$".format(i))
          
          
        
        plt.xlabel(r'$\alpha$')
        plt.ylabel(r'$\sqrt{\lambda}$')
        plt.xlim(0,100)
        plt.legend(loc='upper right')
        plt.grid()
        plt.savefig('wykres1.png')
        
        plt.xscale('log')
        plt.xlim(2,100)
        plt.ylim(0, 0.3)
        plt.savefig('wykres1log.png')
        
    else:
        for i in np.arange(6):
          plt.plot(alfa, vals[:,i], label=r"$u_{{ {0} }}(x)$".format(i))
        plt.xlabel('x')
        plt.ylabel('u(x)')

        plt.legend(loc='upper right')
        plt.grid()
        plt.savefig(f"wykres{filename[-5]}.png")
    


for i in range(1,4):
    plot(f"out{i}.txt")



