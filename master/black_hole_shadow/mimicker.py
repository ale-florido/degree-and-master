import numpy as np
from gr_pyhole import observer, propagator
from gr_pyhole.image import Image
from gr_pyhole.display import Display
from metrica_mimicker import MyMetric

# 1) Set up one of the metrics provided with PyHole
g = MyMetric(1)

# 2) set up an observer
o = observer.Equirectangular(r=15, theta=np.pi/2)

# 3) set up a propagator
p = propagator.SphericalCPU(o, g, Rsky=30)
p.TOLERANCE = 1e-6        # propagator error tolerance

# generate and save an image
i = Image(p,(124,124))   # (128,128) is size in pixels.
i.updateBackground(bgfile='bg-color.png', lines=18)
#i.saveImage('mimi1con5thpiover2.png')
#i.save('mimi1con5thpiover2.npz')
print(str(i))             # prints info about image


d = Display(i, p)         # create new Display
d.show()
