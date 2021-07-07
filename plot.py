#import required package#
import numpy as np
import matplotlib.pyplot as plt

#load#
sal=np.loadtxt('SALPETER.eos',skiprows=1)

#plot#
plt.plot(sal[:,0],sal[:,1],label='Pressure')
plt.plot(sal[:,0],sal[:,2],label='Epsilon',linestyle='--')
plt.plot(sal[:,0],sal[:,3],label='Enthalpy',linestyle='-.')
plt.plot(sal[:,0],sal[:,0],label='Density',linestyle=':')
plt.legend(loc="lower left")
plt.xlabel('Log10 Density (CGS)')
plt.ylabel('Log10 Quantity (CGS)')
plt.title('Plot')
plt.grid(True)
plt.show()
plt.clf()