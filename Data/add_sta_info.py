from glob import glob
from obspy import read
import sys
from obspy.core import AttribDict

argc = len(sys.argv)
if argc != 3:
    print('Usage: python %s sac_dir sta.lst'%sys.argv[0])
    exit(1)

dirname = sys.argv[1]
if dirname[-1] != '/':
    dirname += '/'

sta = {}
with open(sys.argv[2], 'r') as fin:
    for line in fin.readlines():
        tmp = line.strip().split()
        k = tmp[0] + '.' + tmp[1]
        sta[k] = [float(tmp[2]), float(tmp[3])]

sacs = sorted(glob(dirname+'*.SAC'))
for f in sacs:
    print('Add info to', f, '...')
    tr = read(f)[0]
    k = tr.stats.network + '.' + tr.stats.station
    if 'sac' not in tr.stats:
        tr.stats.sac = AttribDict({})
    tr.stats.sac.stlo = sta[k][0]
    tr.stats.sac.stla = sta[k][1]
    tr.stats.sac.lcalda = 1
    tr.write(f, format='SAC')
