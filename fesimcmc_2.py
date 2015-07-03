import pylabsetup

import numpy as np
from scipy.io.idl import readsav
import pylab as pl
from scipy.ndimage import filters
from scipy.signal import  gaussian
from scipy.interpolate import interp1d
from scipy.optimize import minimize 
import time,sys,os
from numpy import size,mean,std,median,where
from pylab import figure

def fittemplate(p,fmean_input,wlog_input,lx,ly,ly_err, plot=False, fig=None):
    #scale2=-p[3]*100
    scale1,v,sig,a=p[2]*100,-p[1]*1000,p[0]*10,p[3]
    w_lower=lx[0]+scale1
    w_upper=lx[-1]-scale1
 
    inds=(lx>w_lower) * (lx<w_upper)
    ly_new=ly[inds]
    ly_err_new=ly_err[inds]
    lx_new=lx[inds]
 
    doppler=1.0+v/299792.458
    b = gaussian(300, sig)
    thisy=filters.convolve1d(a*fmean_input, b/b.sum())  # first do Gaussian convolution
    f2 = interp1d(wlog_input*doppler, thisy,bounds_error=False, fill_value=0)(lx_new) 
    chisq=np.sum((ly_new-f2)**2/ly_err_new**2)/(len(ly_new)-len(p))
    
    if plot:
        if not fig:
            fig=figure(figsize=(10,10))
        fig.add_subplot(1,1,1)
        #pl.plot(wlog_input, fmean_input)
        pl.plot(x_flat,y_flat,'k',alpha=0.5, label="SN PTF10qts")
        pl.plot(lx_new,ly_new,'k')
#        pl.plot(wlog_input*doppler, fmean_input,'r',linewidth=2,alpha=0.5,label=r"scale1=%.2f $\AA$ v=%.2f km/s   $\sigma$=%.2f $\AA$   $\chi^2$=%.2f"%(scale1, -v,sig,chisq))
        pl.plot(wlog_input*doppler, fmean_input,'r',linewidth=2,alpha=0.5,label="template (mean Ic spectrum)")
        pl.plot(lx_new,f2,'r',linewidth=4,label="shifted, broadened template")
        pl.text(5000,0.25,r"v=%.2f km/s   $\sigma$=%.2f $\AA$  $\chi^2$=%.2f"%(-v,sig,chisq),fontsize=18)#%(-v,sig,chisq))
        pl.text(4700,0.18,"Fe",ha='center',fontsize=30)
        pl.text(5700,0.18,"Si",ha='center',fontsize=30)
        pl.text(5000,0.28,"gaussian kernel parameters",fontsize=18)#%(-v,sig,chisq))
#"v=%.2f km/s   $\sigma$=%.2f $\AA$  $\chi^2$=%.2f"%(-v,sig,chisq)
        pl.xlim(3200,8500)
        pl.ylim(-0.35,0.450)
        pl.xlabel("wavelength (A)",fontsize=20)
        pl.ylabel("flux (normalized)",fontsize=20)
        pl.legend(fontsize=15)
    
    #print p,chisq
    return chisq
 

dd=np.loadtxt(sys.argv[1], unpack=True) 
x,y=dd[0],dd[1] # wavelength, flux


s=readsav('meanspecIc_1specperSN_0.sav') # read in Ic template

#pl.plot(s.wlog,s.fmean)
#pl.xlabel("wavelength (A)")
#pl.ylabel("flux (normalized)")
#pl.title ("template spectrum (Ic-mean)")
 

wlog_input=s.wlog[where((s.wlog > 4400) & (s.wlog < 9000))]
fmean_input=s.fmean[where((s.wlog > 4400) & (s.wlog < 9000))]
 
#pl.plot(s.wlog, s.fmean,alpha=0.5) 
#pl.xlim(2000, 10000)
#pl.xlabel("wavelength (A)")
#pl.ylabel("flux (normalized)")
#pl.title ("template spectrum (Ic-mean)")
#pl.figure(1)
#pl.plot(wlog_input, fmean_input,'b',linewidth=3) # plot Ic template - the shift & convolution part
#pl.xlim(2000, 10000)

s2=readsav(sys.argv[1]+"-flat.sav") # read in Icbl flattened spectra
x_flat=s2.wavelog_input[0:1024]
y_flat=s2.flatflux_input
y_flat_err=s2.flatflux_err_input

size(x_flat),size(y_flat),size(y_flat_err), size(y_flat-y_flat_err)

'''pl.figure(figsize=(15,15))
pl.title("PTF 10qts 20100815",fontsize=15)
pl.plot(x,y,'b',label="original spectrum") # plot Icbl raw spectra
pl.xlim(2000, 10000)
pl.plot(x_flat,y_flat,'r',label="flat spectrum") # plot Icbl flattened spectra
pl.xlim(2000, 10000)
#pl.plot(x_flat,y_flat_err,'g',label="errors of flat spectrum") # plot errors of Icbl flattened spectra
#pl.xlim(2000, 10000)
pl.fill_between(x_flat,y_flat-y_flat_err,y_flat+y_flat_err, color='r', alpha=0.2) # plot errors of Icbl flattened spectra
pl.xlabel("wavelength (A)",fontsize=15)
pl.ylabel("flux (normalized)",fontsize=15)
pl.legend() 

pl.figure(figsize=(15,15))
pl.title("PTF 10qts 20100815")
pl.plot(x_flat,y_flat,'k',alpha=0.5) # plot Icbl flattened spectra
'''
Fe_lower=4400 
Fe_upper=5200
Si_lower=5100
Si_upper=6400
Fex=x_flat[np.where((x_flat>Fe_lower) & (x_flat<Fe_upper))]
Fey=y_flat[np.where((x_flat>Fe_lower) & (x_flat<Fe_upper))]
Fey_err=y_flat_err[np.where((x_flat>Fe_lower) & (x_flat<Fe_upper))]
Six=x_flat[np.where((x_flat>Si_lower) & (x_flat<Si_upper))]
Siy=y_flat[np.where((x_flat>Si_lower) & (x_flat<Si_upper))]
Siy_err=y_flat_err[np.where((x_flat>Si_lower) & (x_flat<Si_upper))]
'''pl.plot(Fex,Fey,'k',linewidth=2)
pl.plot(Six,Siy,'k',linewidth=2)
pl.text(mean(Fex),0.2,"Fe",ha='center',fontsize=30)
pl.text(mean(Six),0.2,"Si",ha='center',fontsize=30)
pl.xlim(2000, 10000)
pl.xlabel("wavelength")
pl.ylabel("flux (normalized)")
'''

p0=np.array([2,8,1,1])
#fig=figure(figsize=(5,5))
#fittemplate(p0,fmean_input,wlog_input,Fex,Fey,Fey_err,plot=True,fig=fig)


t1 = time.time()
resFe=minimize(fittemplate, p0, args=(fmean_input,wlog_input,Fex,Fey, Fey_err), method='Powell')
t2 = time.time()
print 'minimization took {} seconds'.format(t2 - t1)
print "Fe minimization with initial guess: sigma=",p0[0]*10,"A, v=", p0[1]*1000,"km/s", p0[2]*100, "A"#, p0[3]*100, "A"
print " sigma:",resFe['x'][0]*10,"A, v=", resFe['x'][1]*1000,"km/s", resFe['x'][2]*100, 'A'#, resFe['x'][3]*100, 'A'
 

p0=np.array([2,3,1,1])
#fig=figure(figsize=(5,5))
#fittemplate(p0,fmean_input,wlog_input,Six,Siy,Siy_err,plot=True,fig=fig)
t1 = time.time()
resSi=minimize(fittemplate, p0, args=(fmean_input,wlog_input,Six,Siy,Siy_err), method='Powell')
t2 = time.time()
print 'minimization took {} seconds'.format(t2 - t1)
print "Si minimization with initial guess: sigma=",p0[0]*10,"A, v=", p0[1]*1000,"km/s", p0[2]*100, "A"#, p0[3]*100, "A"
print " sigma:",resSi['x'][0]*10,"A, v=", resSi['x'][1]*1000,"km/s", resSi['x'][2]*100,"A"#, resSi['x'][3]*100,"A"


p0=resFe['x']
#fig=figure(figsize=(5,5))
#fittemplate(p0,fmean_input,wlog_input,Six,Siy,Siy_err,plot=True,fig=fig)
t1 = time.time()
resSi=minimize(fittemplate, p0, args=(fmean_input,wlog_input,Six,Siy,Siy_err), method='Powell')
t2 = time.time()
print 'minimization took {} seconds'.format(t2 - t1)
print "Si minimization with initial guess: sigma=",p0[0]*10,"A, v=", p0[1]*1000,"km/s", p0[2]*100, "A"#, p0[3]*100, "A"
print " sigma:",resSi['x'][0]*10,"A, v=", resSi['x'][1]*1000,"km/s",resSi['x'][2]*100,"A"#,resSi['x'][3]*100,"A"




p0=np.array([2,8,1,1])
#fig=figure(figsize=(5,5))
#fittemplate(p0,fmean_input,wlog_input,Six,Siy,Siy_err,plot=True,fig=fig)
t1 = time.time()
resSi=minimize(fittemplate, p0, args=(fmean_input,wlog_input,Six,Siy, Siy_err), method='Powell')
t2 = time.time()
print 'minimization took {} seconds'.format(t2 - t1)
print "Si minimization with initial guess: sigma=",p0[0]*10,"A, v=", p0[1]*1000,"km/s", p0[2]*100, "A"#, p0[3]*100, "A"
print " sigma:",resSi['x'][0]*10,"A, v=", resSi['x'][1]*1000,"km/s",resSi['x'][2]*100, "A"#,resSi['x'][3]*100, "A"


p0=resSi['x']
#fig=figure(figsize=(5,5))
#fittemplate(p0,fmean_input,wlog_input,Fex,Fey,Fey_err,plot=True,fig=fig)
t1 = time.time()
resFe=minimize(fittemplate, p0, args=(fmean_input,wlog_input,Fex,Fey, Fey_err), method='Powell')
t2 = time.time()
print 'minimization took {} seconds'.format(t2 - t1)
print "Fe minimization with initial guess: sigma=",p0[0]*10,"A, v=", p0[1]*1000,"km/s",p0[2]*100,"A"#,p0[3]*100,"A"
print " sigma:",resFe['x'][0]*10,"A, v=", resFe['x'][1]*1000,"km/s", "scale: ", resFe['x'][2]*100#, "A", resFe['x'][3]*100, "A"

fittemplate(resFe['x'],fmean_input,wlog_input,Fex,Fey,Fey_err,plot=True)
fittemplate(resSi['x'],fmean_input,wlog_input,Six,Siy,Siy_err,plot=True)

def logprior (p):
    s=p[0]
    v=p[1]
    scale1=p[2]
    #scale2=p[3]
    
    if v>0 and s>0 and v>0.5 and scale1>0 and scale1<2 and p[3]>0: # and scale2>0 and scale1<3 and scale2<3: 
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
##ndim, nwalkers = 2, 32
#fig1=pl.figure(10)
#fig2=pl.figure(11)
ndim, nwalkers = 4, 6*2
 
    
best_pos=[]
#p00=np.array([2,8,1])
#p0 = [p00 + 1e-6*np.random.randn(ndim) for i in range(nwalkers)]
p0 = [resFe["x"] + 1e-6*np.random.randn(ndim) for i in range(nwalkers)]
samplerFe = emcee.EnsembleSampler(nwalkers, ndim, logp, 
                                        args=(Fex,Fey,Fes,fmean_input,wlog_input))
pos, prob, state = samplerFe.run_mcmc(p0, 30)
best_pos.append(samplerFe.flatchain[samplerFe.flatlnprobability.argmax()])
pos = emcee.utils.sample_ball(best_pos[-1], best_pos[-1]/100., size=nwalkers)
samplerFe.reset()
pos, prob, state = samplerFe.run_mcmc(pos, 1000)
best_pos.append(samplerFe.flatchain[samplerFe.flatlnprobability.argmax()])
 
#fig1 = triangle.corner(samplerFe.flatchain, 
#                           truths=best_pos[-1], figure=fig1,labels=["r$\sigma$", "$v$", "$scale1$", "norm"])

print("Mean acceptance fraction: {0:.3f}"
                .format(np.mean(samplerFe.acceptance_fraction)))
print best_pos[-1]

'''
for i in range(ndim):
    pl.figure(figsize=(15,2))
    pl.plot(range(1000),samplerFe.chain[:,:,i].T)
    pl.xlabel("steps")
''' 
samples = samplerFe.chain[:, 50:, :].reshape((-1, ndim))
#samples[:, 2] = np.exp(samples[:, 2])
sigma_mcmc, v_mcmc, scale_mcmc, a_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                             zip(*np.percentile(samples, [16, 50, 84],
                                             axis=0)))
print sigma_mcmc, v_mcmc, scale_mcmc, a_mcmc
fig1 = triangle.corner(samplerFe.flatchain, 
                           truths=[sigma_mcmc[0], v_mcmc[0], scale_mcmc[0], a_mcmc[0]],labels=[r"$\sigma$", "$v$", "$scale1$", "norm"])
 
fittemplate([sigma_mcmc[0], v_mcmc[0], scale_mcmc[0], a_mcmc[0]],fmean_input,wlog_input,Fex,Fey,Fes,plot=True)

best_pos=[]
p0 = [resSi["x"] + 1e-6*np.random.randn(ndim) for i in range(nwalkers)]
samplerSi = emcee.EnsembleSampler(nwalkers, ndim, logp, 
                                        args=(Six,Siy,Sis,fmean_input,wlog_input))
pos, prob, state = samplerSi.run_mcmc(p0, 30)
best_pos.append(samplerSi.flatchain[samplerSi.flatlnprobability.argmax()])
pos = emcee.utils.sample_ball(best_pos[-1], best_pos[-1]/100., size=nwalkers)
samplerSi.reset()
pos, prob, state = samplerSi.run_mcmc(p0, 1000)
best_pos.append(samplerSi.flatchain[samplerSi.flatlnprobability.argmax()])
#tmpb= samplerSi.flatchain.T[0][np.where((samplerSi.flatchain.T[0]<best_pos[-1][0]*10)&(samplerSi.flatchain.T[1]<best_pos[-1][1]*10))]
#tmpv= samplerSi.flatchain.T[1][np.where((samplerSi.flatchain.T[0]<best_pos[-1][0]*10)&(samplerSi.flatchain.T[1]<best_pos[-1][1]*10))]
 
#fig2 = triangle.corner(np.array([tmpb,tmpv]).T, 
#fig2 = triangle.corner(samplerSi.flatchain, 
#                       truths=best_pos[-1], figure=fig2,labels=["r$\sigma$", "$v$", "$scale$", "norm"])
print best_pos[-1]
'''
for i in range(ndim):
    pl.figure(figsize=(15,2))
    pl.plot(range(1000),samplerSi.chain[:,:,i].T)
    #pl.ylim(np.median(samplerSi.chain[:,:,i].T[100])-10,np.median(samplerSi.chain[:,:,i].T[100])+10)
    pl.xlabel("steps")
''' 
#fittemplate(best_pos[-1],fmean_input,wlog_input,Six,Siy,Sis,plot=True)

print("Mean acceptance fraction: {0:.3f}"
                      .format(np.mean(samplerSi.acceptance_fraction)))
        
samples = samplerSi.chain[:, 250:, :].reshape((-1, ndim))
#samples[:, 2] = np.exp(samples[:, 2])
Sisigma_mcmc, Siv_mcmc, Siscale_mcmc, Sia_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                             zip(*np.percentile(samples, [16, 50, 84],
                                             axis=0)))
print Sisigma_mcmc, Siv_mcmc, Siscale_mcmc, Sia_mcmc
fig1 = triangle.corner(samplerSi.flatchain, 
                           truths=[Sisigma_mcmc[0], Siv_mcmc[0], Siscale_mcmc[0], Sia_mcmc[0]],labels=[r"$\sigma$", "$v$", "$scale1$", "norm"])
 
 
fittemplate([Sisigma_mcmc[0], Siv_mcmc[0], Siscale_mcmc[0], Sia_mcmc[0]],fmean_input,wlog_input,Six,Siy,Sis,plot=True)

print "Fe", v_mcmc[1],v_mcmc[0],v_mcmc[2]
print "Si", Siv_mcmc[1],Siv_mcmc[0],Siv_mcmc[2]
pl.show()
