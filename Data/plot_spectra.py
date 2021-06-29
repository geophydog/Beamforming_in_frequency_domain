import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import sys
import matplotlib as mpl

argc = len(sys.argv)
if argc != 3:
    print('Usage: python %s input_file km/deg/m'%sys.argv[0])
    exit(1)

infile = sys.argv[1]
unit = sys.argv[2]
data = []
with open(infile, 'r') as fin:
    for i, line in enumerate(fin.readlines()):
        tmp = line.strip().split()
        if i >= 3 and i <= 8:
            data.append(float(tmp[0]))

PP = np.loadtxt('out.txt', dtype=complex)
PP = abs(PP) ** 2

s = np.linspace(data[2], data[3], len(PP[:, 0]))
if unit =='deg':
    s *= 111.19492664455
elif unit == 'm':
    s /= 1e3
else:
    pass
b = np.linspace(data[4], data[5], len(PP[0])) / 180 * np.pi
baz  = np.linspace(b[0], b[-1], 361)
slow = np.linspace(s[0], s[-1], 201)
inp  = interpolate.interp2d(b, s, PP, kind='cubic')
PP   = abs(inp(baz, slow))
PP /= PP.max()
PP = 10 * np.log10(PP)
col = []
with open('cmap.txt', 'r') as fin:
    for line in fin.readlines():
        col.append(line.strip())
cmap = mpl.colors.LinearSegmentedColormap.from_list('cmap', col, 101)

plt.figure(figsize=(7, 6))
ax1 = plt.subplot(111, projection='polar')
plt.pcolormesh(baz, slow, PP, cmap=cmap, vmin=-10, vmax=0)
cbar = plt.colorbar(shrink=0.4, pad=0.2)
cbar.set_label(r'Beam power (dB)', fontsize=14)
cbar.ax.tick_params(labelsize=15)
ax1.grid(color='gray', ls=(10, (5, 8)), lw=1)
ax1.tick_params(axis='y', colors='gray', labelsize=15)
ax1.set_theta_zero_location('N')
ax1.set_theta_direction(-1)
ax1.plot(b, np.ones(len(b))*s.min(), c='k', lw=2)
ax1.plot(b, np.ones(len(b))*s.max(), c='k', lw=2)
if unit == 'deg':
    ax1.set_xlabel('Slowness (s/deg)', fontsize=20)
elif unit == 'm':
    ax1.set_xlabel('Slowness (s/m)', fontsize=20)
else:
    ax1.set_xlabel('Slowness (s/km)', fontsize=20)
plt.xticks(fontsize=17)
plt.yticks(fontsize=17)
plt.title('%.3f-%.3f Hz'%(data[0], data[1]), fontsize=20, pad=15)
plt.tight_layout() 
plt.show()
