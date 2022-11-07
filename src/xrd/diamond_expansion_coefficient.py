from calendar import c
import numpy as np
import scipy.integrate as integrate
from scipy import interpolate
import matplotlib.pyplot as plt

def forward_expansion_coefficient_model(T):
    Xis = np.array([0.0210,0.3897,3.4447,2.2796],dtype=np.float128)
    Debyes = np.array([225.2,634.0,1365.5,3068.8],dtype=np.float128)
    alpha       =   0
    for Xi, Debye in zip(Xis, Debyes):
        alpha       +=  Xi * (Debye / T )**2 * np.exp(Debye / T) / ((np.exp(Debye / T) - 1)**2) * 1e-6
    return alpha

def lattice_constant(T):
    a0                  = 3.56712
    integral            = integrate.quad(forward_expansion_coefficient_model, 0. ,T,)
    a                   = a0 * np.exp(integral[0])
    return a

if __name__=='__main__':
    fig, ax         = plt.subplots(2)
    T_vals          = np.linspace(289,40000.,1000)
    a               =  np.zeros_like(T_vals)
    for id ,T   in enumerate(T_vals):
        a[id]               = lattice_constant(T)
    #np.savetxt('../../.data/expansion_coefficient_diamond.txt',np.hstack([a[:,np.newaxis],T_vals[:,np.newaxis]]))
    ax[0].plot(T_vals,forward_expansion_coefficient_model(T_vals),label='Expansion coefficient')
    ax2             = ax[0].twinx()
    ax2.plot(T_vals,a,label='Lattice constant',c='r')
    yinterp                       = interpolate.interp1d(a,T_vals)
    lattice                       = np.linspace(3.622587690329867-np.sqrt(0.002445515875843866),3.622587690329867+np.sqrt(0.002445515875843866))
    ax[1].plot(lattice,yinterp(lattice),label='Interpolation',c='g')

    plt.show()