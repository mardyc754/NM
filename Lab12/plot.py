import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline

fig = plt.figure()
data = np.loadtxt("out.txt")

simpson = data[:9]
milne = data[9:]

plt.plot(simpson[3:,0], simpson[3:,1], 'o', label="Simpson")
plt.plot(milne[3:,0], milne[3:,1],'o', label="Milne")
plt.plot(simpson[3:,0], simpson[3:,2], 'o', label="Simpson + ekstrapolacja")
plt.plot(milne[3:,0], milne[3:,2],'o', label="Milne + ekstrapolacja")

simpson_x = milne_x = np.linspace(3,8, 500)
simpson_bspline = make_interp_spline(simpson[3:,0], simpson[3:,1])
milne_bspline = make_interp_spline(milne[3:,0], milne[3:,1])

simpson_bspline_richardson = make_interp_spline(simpson[3:,0], simpson[3:,2])
milne_bspline_richardson = make_interp_spline(milne[3:,0], milne[3:,2])

simpson_y = simpson_bspline(simpson_x)
milne_y = milne_bspline(milne_x)
simpson_y_richardson = simpson_bspline_richardson(simpson_x)
milne_y_richardson = milne_bspline_richardson(milne_x)

plt.plot(simpson_x, simpson_y, "--", color="blue")
plt.plot(milne_x, milne_y, "--", color="orange")
plt.plot(simpson_x, simpson_y_richardson, "--", color="green")
plt.plot(milne_x, milne_y_richardson, "--", color="red")

plt.xlabel("Nr iteracji")
plt.ylabel("Wartość całki")
plt.legend(loc='upper right')
plt.grid()

plt.savefig("Dw0_2.png")
