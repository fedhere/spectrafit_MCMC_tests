# -*- coding: utf-8 -*-
# run a specific spectrum; for Icbl at all phases
import sys
import os
import numpy as np
from scipy.io.idl import readsav
import pylab as pl
from scipy.ndimage import filters
from scipy.signal import  gaussian
from scipy.interpolate import interp1d
from scipy.optimize import minimize 
import time
from matplotlib.backends.backend_pdf import PdfPages
import pylabsetup
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

t1 = time.time()

def fittemplate(p,fmean_input,wlog_input,lx,ly,ly_err, plot=False):
    stre,scale1,v,sig=p[3], p[2]*100,-p[1]*1000,p[0]*10
    w_lower=lx[0]+scale1
    w_upper=lx[-1]-scale1

    inds=(lx>w_lower) * (lx<w_upper)
    ly_new=ly[inds]
    ly_err_new=ly_err[inds]
    lx_new=lx[inds]

    beta=v/299792.458
    doppler=np.sqrt((1+beta)/(1-beta))
    b = gaussian(300, sig)
    thisy=filters.convolve1d(stre*fmean_input, b/b.sum())  # first do Gaussian convolution
    f2 = interp1d(wlog_input*doppler, thisy,bounds_error=False, fill_value=0)(lx_new) 
    chisq=np.sum((ly_new-f2)**2/ly_err_new**2)/(len(ly_new)-len(p))
    
    if plot:
        print "got here"
        fig,ax = pl.subplots()
        minorLocatory   = MultipleLocator(0.02)
        minorLocatorx   = MultipleLocator(100)
        ax.xaxis.set_minor_locator(minorLocatorx)
        ax.yaxis.set_minor_locator(minorLocatory)
        pl.plot(x_flat,y_flat,'k',linewidth=2,label="PTF10qts")
#        pl.plot(wlog_input*doppler, stre*fmean_input,'SteelBlue',linewidth=3,label="SN Ic template (a=%.2f v=%.2f km/s)"%(stre, -v))
        pl.plot(wlog_input*doppler,thisy,'r',linewidth=1)
        pl.plot(lx_new,f2,'r',linewidth=4, label="SN Ic template")
        pl.fill_between(x_flat,y_flat-y_flat_err,y_flat+y_flat_err, color='k', alpha=0.3, label="uncertainty") # plot errors of Icbl flattened spectra
        pl.xlabel("wavelength ($\AA$)")
        pl.ylabel("flattened flux")
        pl.xlim([4000, 5900])
        pl.ylim([-0.32,0.59])
        pl.text(5500,-0.21,"phase 0", fontsize=20)
        pl.text(4100,-0.27,r"(v=%.2f, $km~s^{-1}$, w=%.2f $\AA$, a=%.2f, s=%.2f $\AA$)"%(-v,sig,stre,scale1), fontsize=18)
        pl.legend()
        pl.savefig(sn_spec_name+'_paper_2.pdf')
        
        pl.show()

    return chisq

sn_spec_name='10qts_20100815_Lick_3-m_v1-z.flm-flat.sav'
print sn_spec_name
  
s=readsav('meanspecIc_1specperSN_0.sav') # read in Ic template
wlog_input=s.wlog[np.where((s.wlog > 4400) & (s.wlog < 9000))]
fmean_input=s.fmean[np.where((s.wlog > 4400) & (s.wlog < 9000))]

s2=readsav(sn_spec_name) # read in Icbl spectrum
print s2

x_flat=s2.wavelog_input[0:1024]       
y_flat_sm=s2.flatflux_input_sm
y_flat=s2.flatflux_input
y_flat_err=s2.flatflux_err_input

#pp = PdfPages('/Users/yuqianliu/Desktop/regenerateSNID/Icbl_flat_spec_0/'+sn_spec_name+'_paper_1.pdf')
fig,ax = pl.subplots()
minorLocatory   = MultipleLocator(0.02)
minorLocatorx   = MultipleLocator(100)
ax.xaxis.set_minor_locator(minorLocatorx)
ax.yaxis.set_minor_locator(minorLocatory)
pl.plot(x_flat,y_flat,'k',linewidth=2, label="PTF10qts")
pl.fill_between(x_flat,y_flat-y_flat_err,y_flat+y_flat_err, color='k', alpha=0.3, label="uncertainty") # plot errors of Icbl flattened spectra
pl.plot(wlog_input, fmean_input,'SteelBlue', linewidth=3,label="SN Ic template")
pl.text(5500,-0.21,"phase 0", fontsize=20)
pl.xlabel("wavelength ($\AA$)")
pl.ylabel("relative intensity")
pl.xlim([4000, 5900])
pl.ylim([-0.32,0.59])
pl.legend()
pl.savefig(sn_spec_name+'_paper_1.pdf')

#pl.show()

Fe_lower_inds=(x_flat>4100) * (x_flat<4600)
Fe_lower_x_flat=x_flat[Fe_lower_inds]
Fe_lower_y_flat_sm=y_flat_sm[Fe_lower_inds]
Fe_lower=np.min(Fe_lower_x_flat[np.where(Fe_lower_y_flat_sm == np.max(Fe_lower_y_flat_sm))])-100
Fe_upper_inds=(x_flat>4800) * (x_flat<5200)
Fe_upper_x_flat=x_flat[Fe_upper_inds]
Fe_upper_y_flat_sm=y_flat_sm[Fe_upper_inds]
Fe_upper=np.max(Fe_upper_x_flat[np.where(Fe_upper_y_flat_sm == np.max(Fe_upper_y_flat_sm))])+100

Si_lower_inds=(x_flat>4900) * (x_flat<5300)
Si_lower_x_flat=x_flat[Si_lower_inds]
Si_lower_y_flat_sm=y_flat_sm[Si_lower_inds]
Si_lower=np.min(Si_lower_x_flat[np.where(Si_lower_y_flat_sm == np.max(Si_lower_y_flat_sm))])-100
Si_upper_inds=(x_flat>6200) * (x_flat<6600)
Si_upper_x_flat=x_flat[Si_upper_inds]
Si_upper_y_flat_sm=y_flat_sm[Si_upper_inds]
Si_upper=np.max(Si_upper_x_flat[np.where(Si_upper_y_flat_sm == np.max(Si_upper_y_flat_sm))])+100

Fex=x_flat[np.where((x_flat>Fe_lower) & (x_flat<Fe_upper))]
Fey=y_flat[np.where((x_flat>Fe_lower) & (x_flat<Fe_upper))]
Fey_err=y_flat_err[np.where((x_flat>Fe_lower) & (x_flat<Fe_upper))]
Six=x_flat[np.where((x_flat>Si_lower) & (x_flat<Si_upper))]
Siy=y_flat[np.where((x_flat>Si_lower) & (x_flat<Si_upper))]
Siy_err=y_flat_err[np.where((x_flat>Si_lower) & (x_flat<Si_upper))]

p0Fe=np.array([1,8,1,1])
p0Si=np.array([1,8,1,1])

def logprior (p):
    s=p[0]
    v=p[1]
    scale1=p[2]
    stre=p[3]
    #scale2=p[3]
    if v>0 and s>0 and v>0.5 and scale1>0 and scale1<2.0 and stre>0: # and scale2>0 and scale2<3: 
        return 0.0
    return -np.inf

def logl(p,x,y,s,fmean_input,wlog_input):
    return  -np.log(s)-0.5*(fittemplate(p,fmean_input,wlog_input,x,y,s))
    
def logp(p,x,y,s,fmean_input,wlog_input):
    lgl= logl(p,x,y,s,fmean_input,wlog_input)
    #print lgl
    return np.sum(logprior(p)+lgl)

import triangle 
import emcee

Fes=Fey_err #np.ones(len(Fex))/np.sqrt(len(Fex))
Sis=Siy_err #np.ones(len(Six))/np.sqrt(len(Six))
ndim, nwalkers = 4, 6*2

# Fe MCMC
best_pos=[]
p0 = [p0Fe + 1e-6*np.random.randn(ndim) for i in range(nwalkers)]
samplerFe = emcee.EnsembleSampler(nwalkers, ndim, logp, 
                                        args=(Fex,Fey,Fes,fmean_input,wlog_input))
pos, prob, state = samplerFe.run_mcmc(p0, 30) # run MCMC for 30 steps starting from the tiny ball defined above
best_pos.append(samplerFe.flatchain[samplerFe.flatlnprobability.argmax()])
pos = emcee.utils.sample_ball(best_pos[-1], best_pos[-1]/100., size=nwalkers)
samplerFe.reset()
pos, prob, state = samplerFe.run_mcmc(pos, 1000)
best_pos.append(samplerFe.flatchain[samplerFe.flatlnprobability.argmax()])

f=open(sn_spec_name+'_paper.dat','w')
f.write('for Fe \n')
f.write('initial wavelength range: '+str(Fe_lower)+' '+str(Fe_upper)+'\n')
f.write("Mean acceptance fraction: {0:.3f} \n"
                .format(np.mean(samplerFe.acceptance_fraction)))
print("Mean acceptance fraction: {0:.3f}"
                .format(np.mean(samplerFe.acceptance_fraction)))

y_label=[r"$\sigma/10$", "$v/1000$", "$scale/100$", "$stretch$"]

fittemplate(best_pos[-1],fmean_input,wlog_input,Fex,Fey,Fes,plot=True)

print "here "

sys.exit()

f.write('sigma/10, v/1000, scale/100, stretch \n')
f.write('initial guess: '+str(p0Fe)+'\n')
f.write('best value: '+str(best_pos[-1])+'\n')
f.write('16th, 50th, 84th percentiles \n')
f.write(str(np.percentile(samplerFe.chain[:,:,0],[16,50,84]))+'\n')
f.write(str(np.percentile(samplerFe.chain[:,:,1],[16,50,84]))+'\n')
f.write(str(np.percentile(samplerFe.chain[:,:,2],[16,50,84]))+'\n')
f.write(str(np.percentile(samplerFe.chain[:,:,3],[16,50,84]))+'\n')

print best_pos[-1] # sigma/10, velocity/1000, scale/100, stretch
print np.percentile(samplerFe.chain[:,:,0],[16,50,84]) # 16th, 50th, 84th percentiles of the sigma/10 => 1 sigma error bar
print np.percentile(samplerFe.chain[:,:,1],[16,50,84]) # 16th, 50th, 84th percentiles of the velocity/1000
print np.percentile(samplerFe.chain[:,:,2],[16,50,84]) # 16th, 50th, 84th percentiles of the scale/100
print np.percentile(samplerFe.chain[:,:,3],[16,50,84]) # 16th, 50th, 84th percentiles of the stretch

# Si MCMC
best_pos=[]
p0 = [p0Si + 1e-6*np.random.randn(ndim) for i in range(nwalkers)]
samplerSi = emcee.EnsembleSampler(nwalkers, ndim, logp, 
                                        args=(Six,Siy,Sis,fmean_input,wlog_input))
pos, prob, state = samplerSi.run_mcmc(p0, 30)
best_pos.append(samplerSi.flatchain[samplerSi.flatlnprobability.argmax()])
pos = emcee.utils.sample_ball(best_pos[-1], best_pos[-1]/100., size=nwalkers)
samplerSi.reset()
pos, prob, state = samplerSi.run_mcmc(p0, 1000)
best_pos.append(samplerSi.flatchain[samplerSi.flatlnprobability.argmax()])
    
f.write('for Si \n')
f.write('initial wavelength range: '+str(Si_lower)+' '+str(Si_upper)+'\n')
f.write("Mean acceptance fraction: {0:.3f} \n"
                .format(np.mean(samplerSi.acceptance_fraction)))
print("Mean acceptance fraction: {0:.3f}"
                .format(np.mean(samplerSi.acceptance_fraction)))

fittemplate(best_pos[-1],fmean_input,wlog_input,Six,Siy,Sis,plot=True)
    
f.write('sigma/10, v/1000, scale/100, stretch \n')
f.write('initial guess: '+str(p0Si)+'\n')   
f.write('best value: '+str(best_pos[-1])+'\n')
f.write('16th, 50th, 84th percentiles \n')
f.write(str(np.percentile(samplerSi.chain[:,:,0],[16,50,84]))+'\n')
f.write(str(np.percentile(samplerSi.chain[:,:,1],[16,50,84]))+'\n')
f.write(str(np.percentile(samplerSi.chain[:,:,2],[16,50,84]))+'\n')
f.write(str(np.percentile(samplerSi.chain[:,:,3],[16,50,84]))+'\n')

f.close()

print best_pos[-1] # sigma/10, velocity/1000, scale/100 
print np.percentile(samplerSi.chain[:,:,0],[16,50,84]) # 16th, 50th, 84th percentiles of the sigma/10
print np.percentile(samplerSi.chain[:,:,1],[16,50,84]) # 16th, 50th, 84th percentiles of the velocity/1000
print np.percentile(samplerSi.chain[:,:,2],[16,50,84]) # 16th, 50th, 84th percentiles of the scale/100
print np.percentile(samplerSi.chain[:,:,3],[16,50,84]) # 16th, 50th, 84th percentiles of the stretch

t2 = time.time()

print ' '
print p0Fe
print Fe_lower, Fe_upper
print p0Si
print Si_lower, Si_upper
print 'minimization took {} seconds'.format(t2 - t1)

pp.close()

