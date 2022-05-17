import numpy as np
import matplotlib.pyplot as plt
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

p = np.loadtxt('out.txt', dtype=complex)
p = np.abs(p)
p /= p.max()
p = 10 * np.log10(p**2)
ns, nb = PP.shape
b = np.linspace(data[4], data[5], nb) / 180 * np.pi
s = np.linspace(data[2], data[3], ns)

if unit =='deg':
    s *= 111.19492664455
elif unit == 'm':
    s /= 1e3
else:
    pass

with open('cmap.txt', 'r') as fin:
    col = fin.readline().split()
cmap = mpl.colors.LinearSegmentedColormap.from_list('cmap', col, 201)

plt.figure(figsize=(7, 6))
ax1 = plt.subplot(111, projection='polar')
plt.pcolormesh(b, s, p, cmap=cmap, vmin=-10, vmax=0)
cbar = plt.colorbar(shrink=0.4, pad=0.2)
cbar.set_label(r'Beam power (dB)', fontsize=14)
cbar.ax.tick_params(labelsize=15)
ax1.grid(color='gray', ls=(10, (5, 8)), lw=1)
ax1.tick_params(axis='y', colors='gray', labelsize=15)
ax1.set_theta_zero_location('N')
ax1.set_theta_direction(-1)
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
