import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure()

file = "out.txt"
data = np.loadtxt(file)

t = data[:, 0]	
x = data[:, 1]

plt.plot(t, x, 'bo', lw=2, label='x(t), h=0.1')

l1 = plt.legend()
plt.grid()

plt.xlabel('t')
plt.ylabel('x(t)')
plt.title('wychylenie x(t)')

plt.savefig("zad1.png")

