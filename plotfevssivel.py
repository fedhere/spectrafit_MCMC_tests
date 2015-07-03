import os,inspect,sys
cmd_folder = os.getenv("SESNCFAlib")
import numpy as np

if cmd_folder not in sys.path:
         sys.path.insert(0, cmd_folder)
    

from pylabsetup import *                                                                                                     

f=open("velocitiesFevsSi.dat")

f.readline()
#'SN name\tSi_shift (km/s)\tSi_new (km/s)\tFe_shift (km/s)\tFe_new (km/s)\n'

sn,si,fe,sim,sip,fem,fep=[],[],[],[],[],[],[]

for l in f:
    if l.startswith('#'):
         continue
    print l
    l=l.strip().split()
    sn.append(l[0])
    fem.append(float(l[1]))
    fe.append(float(l[2]))
    fep.append(float(l[3]))
    sim.append(float(l[4]))
    si.append(float(l[5]))
    sip.append(float(l[6]))
    


fe=np.array(fe)
si=np.array(si)

sierr=[si-np.array(sim),np.array(sip)-si]
feerr=[fe-np.array(fem),np.array(fep)-fe]
print sierr
print feerr
sys.exit()
#pl.figure(figsize=(15,15))
ax=pl.subplot(121, aspect='equal')
ax.minorticks_on()
#pl.ylim(9,29)
#pl.xlim(9,29)
pl.errorbar(si,fe,xerr=sierr,yerr=feerr,fmt= 'ko')
pl.plot([min(si),max(fe)],[min(si),max(fe)], 'k-')
pl.xlabel ("Si II 6355 (1000 km/s)")
pl.ylabel ("Fe II 5169 (1000 km/s)")
pl.title ("Ic-BL Photospheric velocity")
pl.savefig("fevssivelocities.png",bbox_inches='tight', dpi=150)

pl.show()
