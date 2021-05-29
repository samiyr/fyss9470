from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
pp   = PdfPages('azimuth_angles.pdf')
tmp1 = plt.figure(1)
tmp1.set_size_inches(8.00,6.00)
plot = open('azimuth_angles-0.dat')
plot = [line.split() for line in plot]
valx = [float(x[0]) for x in plot]
valy = [float(x[1]) for x in plot]
vale = [float(x[2]) for x in plot]
plt.hist( valx, vale, weights = valy, histtype='step', label=r'azimuth angles')
plt.xlim( -4.000e+00, 4.000e+00)
plt.ylim( 0.000e+00, 1.312e+02)
plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,3))
plt.legend(frameon=False,loc='best')
plt.title(r'Azimuth angles')
plt.xlabel(r'$\phi$')
plt.ylabel(r'N')
pp.savefig(tmp1,bbox_inches='tight')
plt.clf()
pp.close()
