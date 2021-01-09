import numpy as np

n = 60
r = 3
dist = 50

print('%2d              # Number of stations (nsta)'%n)
for i in range(n):
    xy = np.random.randn(2)
    x = xy[0] * r + dist
    y = xy[1] * r
    print('%4.1f  %4.1f   STA%02d       # Station Lat Lon Name'%(y, x, i+1))

