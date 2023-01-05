import numpy as np
import pylab as plt

ref = np.loadtxt('res.txt',dtype="float64")
y = ref[:,2]
l = int(np.sqrt(len(y)))
y_grid = y.reshape((l,l))

plt.imshow(y_grid, cmap="hot")
plt.gca().invert_yaxis()
plt.axis(False)
plt.colorbar()

plt.show()