import numpy as np
from scipy.io.idl import readsav
import pylab as pl
from scipy.ndimage import filters
from scipy.signal import  gaussian
from scipy.interpolate import interp1d
from scipy.optimize import minimize 
import time
from numpy import mean,median,std,size
#%pylab inline
pl.ion()
def fittemplate(p,fmean_input,wlog_input,lx,ly,ly_err, plot=False):
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
        pl.figure(figsize=(15,15))
        #pl.plot(wlog_input, fmean_input)
        pl.plot(x_flat,y_flat,'k',alpha=0.5)
        pl.plot(lx_new,ly_new,'k')
        pl.plot(lx_new,f2,'r',linewidth=3)
        pl.plot(wlog_input*doppler, fmean_input,'r',linewidth=2,alpha=0.5,label=r"scale1=%.2f $\AA$ v=%.2f km/s   $\sigma$=%.2f $\AA$   $\chi^2$=%.2f"%(scale1, -v,sig,chisq))
        pl.xlabel("wavelength (A)",fontsize=15)
        pl.ylabel("flux (normalized)",fontsize=15)
        pl.legend(fontsize=15)
        pl.show()
        pl.figure(figsize=(5,5))
        #pl.plot(wlog_input, fmean_input)
        pl.xlim(3200,8400)
        pl.ylim(-0.75,0.62)
        pl.plot(x_flat,y_flat,'k',alpha=0.5, label="SN2002ap, phase 0")
        pl.plot(lx_new,ly_new,'k')
        pl.plot(lx_new,f2,'r',linewidth=3, label="raw Ic template")
        pl.plot(wlog_input*doppler, fmean_input,'r',linewidth=2,alpha=0.5,label=r"convolved Ic template: v=%.2f km/s   $\sigma$=%.2f $\AA$"%(-v,sig))
        pl.xlabel("wavelength (A)",fontsize=18)
        pl.ylabel("flux (normalized)",fontsize=18)
        pl.legend(fontsize=18)
        pl.draw()
        raw_input()
        pl.show()
    #print p,chisq
    return chisq
dd=np.loadtxt("10qts_20100815_Lick_3-m_v1-z.flm", unpack=True) 
x,y=dd[0],dd[1] # wavelength, flux
s=readsav('meanspecIc_1specperSN_0.sav') # read in Ic template
pl.plot(s.wlog,s.fmean)
pl.xlabel("wavelength (A)")
pl.ylabel("flux (normalized)")
pl.title ("template spectrum (Ic-mean)")
pl.show()
wlog_input=s.wlog[(s.wlog > 4400) * (s.wlog < 9000)]
fmean_input=s.fmean[(s.wlog > 4400) * (s.wlog < 9000)]

pl.plot(s.wlog, s.fmean,alpha=0.5) 
pl.xlim(2000, 10000)
pl.xlabel("wavelength (A)")
pl.ylabel("flux (normalized)")
pl.title ("template spectrum (Ic-mean)")
pl.figure(1)
pl.plot(wlog_input, fmean_input,'b',linewidth=3) # plot Ic template - the shift & convolution part
pl.xlim(2000, 10000)
pl.show()

#s2=readsav('10qts_20100815_Lick_3-m_v1-z.flm-flat.sav') # read in Icbl flattened spectra
s2=readsav('sn2002ap-20020208-z.flm-flat.sav')
#s2=readsav('sn2003jd-20031101-z-bl.flm-flat.sav')
#s2=readsav('sn2007ru-20071206-z.flm-flat.sav')
x_flat=s2.wavelog_input[0:1024]
y_flat=s2.flatflux_input
y_flat_err=s2.flatflux_err_input
size(x_flat),size(y_flat),size(y_flat_err), size(y_flat-y_flat_err)

pl.figure(figsize=(15,15))
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
pl.show()
pl.figure(figsize=(15,15))
pl.title("PTF 10qts 20100815")
pl.plot(x_flat,y_flat,'k',alpha=0.5) # plot Icbl flattened spectra
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
pl.plot(Fex,Fey,'k',linewidth=2)
pl.plot(Six,Siy,'k',linewidth=2)
pl.text(mean(Fex),0.2,"Fe",ha='center',fontsize=30)
pl.text(mean(Six),0.2,"Si",ha='center',fontsize=30)
pl.xlim(2000, 10000)
pl.xlabel("wavelength")
pl.ylabel("flux (normalized)")
pl.show()

p0=np.array([2,8,1,1])
fittemplate(p0,fmean_input,wlog_input,Fex,Fey,Fey_err)

t1 = time.time()
resFe=minimize(fittemplate, p0, args=(fmean_input,wlog_input,Fex,Fey, Fey_err), method='Powell')
t2 = time.time()
print 'minimization took {} seconds'.format(t2 - t1)
print "Fe minimization with initial guess: sigma=",p0[0]*10,"A, v=", p0[1]*1000,"km/s", p0[2]*100, "A"#, p0[3]*100, "A"
print " sigma:",resFe['x'][0]*10,"A, v=", resFe['x'][1]*1000,"km/s", resFe['x'][2]*100, 'A'#, resFe['x'][3]*100, 'A'

p0=np.array([2,3,1,1])
fittemplate(p0,fmean_input,wlog_input,Six,Siy,Siy_err)
t1 = time.time()
resSi=minimize(fittemplate, p0, args=(fmean_input,wlog_input,Six,Siy,Siy_err), method='Powell')
t2 = time.time()
print 'minimization took {} seconds'.format(t2 - t1)
print "Si minimization with initial guess: sigma=",p0[0]*10,"A, v=", p0[1]*1000,"km/s", p0[2]*100, "A"#, p0[3]*100, "A"
print " sigma:",resSi['x'][0]*10,"A, v=", resSi['x'][1]*1000,"km/s", resSi['x'][2]*100,"A"#, resSi['x'][3]*100,"A"

p0=resFe['x']
fittemplate(p0,fmean_input,wlog_input,Six,Siy,Siy_err)
t1 = time.time()
resSi=minimize(fittemplate, p0, args=(fmean_input,wlog_input,Six,Siy,Siy_err), method='Powell')
t2 = time.time()
print 'minimization took {} seconds'.format(t2 - t1)
print "Si minimization with initial guess: sigma=",p0[0]*10,"A, v=", p0[1]*1000,"km/s", p0[2]*100, "A"#, p0[3]*100, "A"
print " sigma:",resSi['x'][0]*10,"A, v=", resSi['x'][1]*1000,"km/s",resSi['x'][2]*100,"A"#,resSi['x'][3]*100,"A"

p0=np.array([2,8,1,1])
fittemplate(p0,fmean_input,wlog_input,Six,Siy,Siy_err)
t1 = time.time()
resSi=minimize(fittemplate, p0, args=(fmean_input,wlog_input,Six,Siy, Siy_err), method='Powell')
t2 = time.time()
print 'minimization took {} seconds'.format(t2 - t1)
print "Si minimization with initial guess: sigma=",p0[0]*10,"A, v=", p0[1]*1000,"km/s", p0[2]*100, "A"#, p0[3]*100, "A"
print " sigma:",resSi['x'][0]*10,"A, v=", resSi['x'][1]*1000,"km/s",resSi['x'][2]*100, "A"#,resSi['x'][3]*100, "A"

p0=resSi['x']
fittemplate(p0,fmean_input,wlog_input,Fex,Fey,Fey_err)
t1 = time.time()
resFe=minimize(fittemplate, p0, args=(fmean_input,wlog_input,Fex,Fey, Fey_err), method='Powell')
t2 = time.time()
print 'minimization took {} seconds'.format(t2 - t1)
print "Fe minimization with initial guess: sigma=",p0[0]*10,"A, v=", p0[1]*1000,"km/s",p0[2]*100,"A"#,p0[3]*100,"A"
print " sigma:",resFe['x'][0]*10,"A, v=", resFe['x'][1]*1000,"km/s", "scale: ", resFe['x'][2]*100#, "A", resFe['x'][3]*100, "A"

fittemplate(resFe['x'],fmean_input,wlog_input,Fex,Fey,Fey_err,plot=True)
fittemplate(resSi['x'],fmean_input,wlog_input,Six,Siy,Siy_err,plot=True)

fittemplate([1.74,8.430,1,1],fmean_input,wlog_input,Fex,Fey,Fey_err,plot=True)
fittemplate([1.09,4.919,1,1],fmean_input,wlog_input,Six,Siy,Siy_err,plot=True)

def logprior (p):
    s=p[0]
    v=p[1]
    scale1=p[2]
    #scale2=p[3]
    if v>0 and s>0 and v>0.5 and scale1>0 and scale1<3: # and scale2>0 and scale1<3 and scale2<3: 
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
fig1=pl.figure(10)
fig2=pl.figure(11)
ndim, nwalkers = 3, 6*2

    
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

fig1 = triangle.corner(samplerFe.flatchain, 
                           truths=best_pos[-1], figure=fig1,labels=["r$\sigma$", "$v$", "$scale1$"])

print("Mean acceptance fraction: {0:.3f}"
                .format(np.mean(samplerFe.acceptance_fraction)))
print best_pos[-1]
for i in range(ndim):
    pl.figure(figsize=(15,2))
    pl.plot(range(1000),samplerFe.chain[:,:,i].T)
    pl.xlabel("steps")

fittemplate(best_pos[-1],fmean_input,wlog_input,Fex,Fey,Fes,plot=True)

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
fig2 = triangle.corner(samplerSi.flatchain, 
                       truths=best_pos[-1], figure=fig2,labels=["r$\sigma$", "$v$", "$scale$"])
print("Mean acceptance fraction: {0:.3f}"
                .format(np.mean(samplerSi.acceptance_fraction)))
print best_pos[-1]
for i in range(ndim):
    pl.figure(figsize=(15,2))
    pl.plot(range(1000),samplerSi.chain[:,:,i].T)
    #pl.ylim(np.median(samplerSi.chain[:,:,i].T[100])-10,np.median(samplerSi.chain[:,:,i].T[100])+10)
    pl.xlabel("steps")

fittemplate(best_pos[-1],fmean_input,wlog_input,Six,Siy,Sis,plot=True)
