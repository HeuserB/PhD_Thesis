import numpy as np

def gaussian(x, amp_gauss, cen_gauss, wid_gauss):
#    "1-d gaussian: gaussian(x, amp, cen, wid)"
      return (amp_gauss/(np.sqrt(2*np.pi)*wid_gauss)) * np.exp(-(x-cen_gauss)**2 /(2*wid_gauss**2))

def background(x,a,b):
        return a * x + b

def lingaussian(x,slope,intercept,amp_gauss, cen_gauss, wid_gauss):
#    "1-d gaussian: lingaussian(x, slope,intercept,amp, cen, wid)"
    return (slope*x+intercept+(amp_gauss/(np.sqrt(2*np.pi)*wid_gauss)) * np.exp(-(x-cen_gauss)**2 /(2*wid_gauss**2)))

def lin_gauss_lorentzian(x, slope,intercept,amp_gauss, cen_gauss, wid_gauss,amp_lor, cen_lor, wid_lor):
#    "1-d lorentzian: lorentzian(x, amp, cen, wid)"
    return ((slope*x+intercept)\
              +(amp_gauss/(np.sqrt(2*np.pi)*wid_gauss)) * np.exp(-(x-cen_gauss)**2 /(2*wid_gauss**2))\
              +(amp_lor/np.pi * (wid_lor/((x-cen_lor)**2+wid_lor**2))))

def lin_gauss_2lorentzian(x, slope,intercept,amp_gauss, cen_gauss, wid_gauss, amp_lor, cen_lor, wid_lor, amp_lor2, cen_lor2, wid_lor2):
#    "1-d lorentzian: lorentzian(x, amp, cen, wid)"
    return ((slope*x+intercept)\
              +(amp_gauss/(np.sqrt(2*np.pi)*wid_gauss)) * np.exp(-(x-cen_gauss)**2 /(2*wid_gauss**2))\
              +(amp_lor/np.pi * (wid_lor/((x-cen_lor)**2+wid_lor**2))))\
              +(amp_lor2/np.pi * (wid_lor2/((x-cen_lor2)**2+wid_lor2**2)))

def lorentzian(x,amp_lor,cen_lor,wid_lor):
    return (amp_lor/np.pi * (wid_lor/((x-cen_lor)**2+wid_lor**2)))

def gauss_lorentzian(x,amp_gauss, cen_gauss, wid_gauss,amp_lor,cen_lor,wid_lor):
    return (amp_lor/np.pi * (wid_lor/((x-cen_lor)**2+wid_lor**2))) + (amp_gauss/(np.sqrt(2*np.pi)*wid_gauss)) * np.exp(-(x-cen_gauss)**2 /(2*wid_gauss**2))