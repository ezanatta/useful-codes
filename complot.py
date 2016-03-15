import numpy as np
import pylab

r = np.loadtxt('PNecompletness.dat', usecols=(0,))
c = np.loadtxt('PNecompletness.dat', usecols=(1,))


fit = np.polyfit(r, c, 6)
p = np.poly1d(fit)

pol = []

for i in range(0,400):
    pol.append(p(i))

x = np.linspace(0, np.amax(r),400)

rad = np.loadtxt('PNEdensity-circularbins.dat', usecols=(6,))
prad = []
for i in rad:
    prad.append(p(i))

with open('PNeextrapolation.dat', 'w') as w:
    for i in range(0,len(rad)):
        print >>w, rad[i], prad[i]

pylab.plot(r, c, 'bo', label='Lodo Completeness')
pylab.xlabel('R (arcseconds)')
pylab.ylabel('% of Completness')
pylab.plot(x, pol, label='Extrapolation')
pylab.legend(loc='lower right')
pylab.show()