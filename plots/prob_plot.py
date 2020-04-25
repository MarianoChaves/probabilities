import numpy as np
import matplotlib.pyplot as plt
import math

En,Numerical,My,SO,MySO=np.loadtxt("prob.dat",dtype=float,delimiter="\t",usecols=(0,1,2,3,4),comments="#",unpack=True)



fig = plt.figure()
axes=fig.add_subplot(1,1,1, axisbg="1.0")


axes.plot(En,SO,color='cyan',label='Standard Oscillation')
axes.plot(En,Numerical,color='black',label='Numerical')
axes.scatter(En,My,color='black',label='This Work',s=5)

axes.legend()

y_axis_NAME='$P(\\nu_\\mu\\rightarrow\\nu_e)$'
x_axis_NAME='$E_\\nu(GeV)$'

axes.set_xscale('log')
axes.set_xlim([0.4,10.0])
axes.set_ylim([0.0,0.1])
axes.set_ylabel(y_axis_NAME)
axes.set_xlabel(x_axis_NAME)

plt.title('DUNE case, $\\rho = 2.8 $g/cm$^3$ and 1300 km')
plt.savefig('prob.pdf') 
plt.show()




