import numpy as np
import matplotlib.pyplot as plt

def fitPlaneLTSQ(XYZ):
    A = np.c_[XYZ[:,0], XYZ[:,1], np.ones(XYZ.shape[0])]
    C,_,_,_ = np.linalg.lstsq(A, XYZ[:,2])    # coefficients
    return C


from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

#image       = cv2.imread('/home/benjamin/Nextcloud/Work/W_PhD/W_PhD_Analysis/AL_samples/Al_probe_0.05/Al_probe_0.05_Nr01.tif')
#points      = np.asarray(image,dtype = np.float64)

points      = np.loadtxt('/home/benjamin/Nextcloud/Work/W_PhD/W_PhD_Analysis/AL_samples/Al_probe_0.05/Al_probe_0.05_Nr01.xyz')
#points      = points.T

print(points.shape)

data = np.c_[points[:,0],points[:,1],points[:,2]]

C           = fitPlaneLTSQ(data)

X,Y         = np.meshgrid(np.linspace(points[:,0].min(),points[:,0].max(),20), \
                np.linspace(points[:,1].min(),points[:,1].max(),20) )

Z = C[0]*X + C[1]*Y + C[2]
print(Z.shape)
points[:,2] -= C[0]*points[:,0] + C[1]*points[:,1] + C[2]

fig         = plt.figure()
ax          = fig.gca(projection='3d')
ax.scatter(points[::50,0], points[::50,1], points[::50,2], s=0.11)
#ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.2)

plt.show()


'''

# Window name in which image is displayed
window_name = 'image'
  
# Using cv2.imshow() method 
# Displaying the image 
cv2.imshow(window_name, image)
  
#waits for user to press any key 
#(this is necessary to avoid Python kernel form crashing)
cv2.waitKey(0) 
  
#closing all open windows 
cv2.destroyAllWindows()
'''