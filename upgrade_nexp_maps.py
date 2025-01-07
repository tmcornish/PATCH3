import healsparse as hsp
import glob

for fd in ['hectomap', 'spring', 'autumn']:
	PATH_MAPS = f'out/{fd}/systmaps/'
	for f in glob.glob(PATH_MAPS + 'decasu_1024_?_nexp_sum.hsp'):
		m = hsp.HealSparseMap.read(f).upgrade(2048)
		m.write(f.replace('1024', '2048'))