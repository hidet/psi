import sys
import time
import datetime
import numpy as np
import pylab as plt
import mass
import psidata
import pyroot as a
from matplotlib.backends.backend_pdf import PdfPages
#import exafs

import cPickle
from os import path

dir_d = "/Users/tatsuno/work/heates/data"
date_p, ds_p = "20141105","O"
date_n, ds_n = "20141105","S"
dir_pl = "%s_%s"%(date_p,ds_p)
ana = "Mn"

# --- booleans ---
forceNew = False  # update
anado    = False  # do analaysis
newcal   = False  # update calibration
resol    = False  # print resolutions
plot     = False  # plot figures for a ch
dump     = False  # dump to ROOT file
# ----------------


# --- get noise date, name, and analysis line
date_n, ds_n, ana = psidata.getnoiseana(dir_pl)        
if ana=="":
    print "Error: ana is empty"
    sys.exit(0)
dir_p = "%s/%s/ljh_files/%s_%s"%(dir_d,date_p,date_p,ds_p)
dir_n = "%s/%s/ljh_files/%s_%s"%(dir_d,date_n,date_n,ds_n)
print "----- pulse: %s "%(dir_pl)
print "----- noise: %s_%s"%(date_n,ds_n)
print "----- ana:   %s"%(ana)
print "----- (forceNew,anado)=(%r,%r)"%(forceNew,anado)
print "----- (newcal,plot,dump)=(%r,%r,%r)"%(newcal,plot,dump)

# --- define basic cuts
pave_high=20000.
if ana=="Se" or ana=="Ge" or ana=="GaAs" or ana=="RbBr" or ana=="highE": pave_high=40000.
cuts = mass.core.controller.AnalysisControl(
    pulse_average=(0, pave_high),
    pretrigger_rms=(1, 100),
    pretrigger_mean_departure_from_median=(-120, 120),
    peak_value=(0, None),
    postpeak_deriv=(None, 250),
    rise_time_ms=(0.10, 0.30),
    peak_time_ms=(0.10, 0.50)
)

# --- set channel numbers
#chans = [i for i in xrange(1,480,2)]
chans = [181,307]
exceptch = [83,109,111,115,185,243,285,323,339,349,351,363,397,401,421]
badch = [245,247,275,277]
exceptch.extend(badch)
for i in exceptch:
    if i in chans: chans.remove(i)
if len(chans)==0:
    raise ValueError("zero number of channels selected")
    sys.exit(0)

# --- load data
data1 = a.load_psi(date_p,"O",date_n,ds_n,chs=chans,cut=cuts,ddir=dir_d)
data2 = a.load_psi(date_p,"P",date_n,ds_n,chs=chans,cut=cuts,ddir=dir_d)
data3 = a.load_psi(date_p,"Q",date_n,ds_n,chs=chans,cut=cuts,ddir=dir_d)
data4 = a.load_psi(date_p,"R",date_n,ds_n,chs=chans,cut=cuts,ddir=dir_d)

dir_p1 = "%s/%s/ljh_files/%s_%s"%(dir_d,date_p,date_p,"O")
dir_p2 = "%s/%s/ljh_files/%s_%s"%(dir_d,date_p,date_p,"P")
dir_p3 = "%s/%s/ljh_files/%s_%s"%(dir_d,date_p,date_p,"Q")
dir_p4 = "%s/%s/ljh_files/%s_%s"%(dir_d,date_p,date_p,"R")
dir_pl1 = "%s_%s"%(date_p,"O")
dir_pl2 = "%s_%s"%(date_p,"P")
dir_pl3 = "%s_%s"%(date_p,"Q")
dir_pl4 = "%s_%s"%(date_p,"R")

# --------------------------------------------------
params = {'xtick.labelsize': 10, # x ticks
          'ytick.labelsize': 10, # y ticks
          'legend.fontsize': 7
         }
plt.rcParams.update(params)

ds1 = data1.datasets[0];  g1 = ds1.good();
ds2 = data2.datasets[0];  g2 = ds2.good();
ds3 = data3.datasets[0];  g3 = ds3.good();
ds4 = data4.datasets[0];  g4 = ds4.good();

nsamples=ds1.nSamples
presamples=ds1.nPresamples
timebase=ds1.timebase
dt = (np.arange(nsamples)-presamples)*timebase*1e3
lt=240
ht=310


import h5py
import commands
# --- read HDF5 file ---
h1 = h5py.File("%s/%s_mass.hdf5"%(dir_p1,dir_pl1))
h2 = h5py.File("%s/%s_mass.hdf5"%(dir_p2,dir_pl2))
h3 = h5py.File("%s/%s_mass.hdf5"%(dir_p3,dir_pl3))
h4 = h5py.File("%s/%s_mass.hdf5"%(dir_p4,dir_pl4))

# ----------------------------------------
ch1 = h1['chan%d'%ds1.channum]
filt1 = ch1['filters']['filt_noconst']
rms1 = filt1.attrs['variance']**0.5
avepulse1 = ch1['average_pulse']
noisepsd1 = ch1['noise_psd']
df1 = noisepsd1.attrs['delta_f']
premean1 = ch1['pretrig_mean']
prerms1 = ch1['pretrig_rms']
timestamp1 = ch1['timestamp']
pulseave1 = ch1['pulse_average']
pulserms1 = ch1['pulse_rms']
filtv1 = ch1['filt_value']
filtvdc1 = ch1['filt_value_dc']
filtvphc1 = ch1['filt_value_phc']
filtvtdc1 = ch1['filt_value_tdc']

ch2 = h2['chan%d'%ds2.channum]
filt2 = ch1['filters']['filt_noconst']
rms2 = filt2.attrs['variance']**0.5
avepulse2 = ch2['average_pulse']
noisepsd2 = ch2['noise_psd']
df2 = noisepsd2.attrs['delta_f']
premean2 = ch2['pretrig_mean']
prerms2 = ch2['pretrig_rms']
timestamp2 = ch2['timestamp']
pulseave2 = ch2['pulse_average']
pulserms2 = ch2['pulse_rms']
filtv2 = ch2['filt_value']
filtvdc2 = ch2['filt_value_dc']
filtvphc2 = ch2['filt_value_phc']
filtvtdc2 = ch2['filt_value_tdc']

ch3 = h3['chan%d'%ds3.channum]
filt3 = ch1['filters']['filt_noconst']
rms3 = filt3.attrs['variance']**0.5
avepulse3 = ch3['average_pulse']
noisepsd3 = ch3['noise_psd']
df3 = noisepsd3.attrs['delta_f']
premean3 = ch3['pretrig_mean']
prerms3 = ch3['pretrig_rms']
timestamp3 = ch3['timestamp']
pulseave3 = ch3['pulse_average']
pulserms3 = ch3['pulse_rms']
filtv3 = ch3['filt_value']
filtvdc3 = ch3['filt_value_dc']
filtvphc3 = ch3['filt_value_phc']
filtvtdc3 = ch3['filt_value_tdc']


ch4 = h4['chan%d'%ds4.channum]
filt4 = ch1['filters']['filt_noconst']
rms4 = filt4.attrs['variance']**0.5
avepulse4 = ch4['average_pulse']
noisepsd4 = ch4['noise_psd']
df4 = noisepsd4.attrs['delta_f']
premean4 = ch4['pretrig_mean']
prerms4 = ch4['pretrig_rms']
timestamp4 = ch4['timestamp']
pulseave4 = ch4['pulse_average']
pulserms4 = ch4['pulse_rms']
filtv4 = ch4['filt_value']
filtvdc4 = ch4['filt_value_dc']
filtvphc4 = ch4['filt_value_phc']
filtvtdc4 = ch4['filt_value_tdc']


# ----------------------------------------

#plt.figure(figsize=(11.69,8.27))
plt.figure(figsize=(8.27,11.69))

# --- Average Pulse ---
ax = plt.subplot(4,2,1)
plt.errorbar(dt, avepulse1, fmt='-', label="Data %s"%dir_pl1, ms=1)                  
plt.errorbar(dt, avepulse2, fmt='-', label="Data %s"%dir_pl2, ms=1)
plt.errorbar(dt, avepulse3, fmt='-', label="Data %s"%dir_pl3, ms=1)
plt.errorbar(dt, avepulse4, fmt='-', label="Data %s"%dir_pl4, ms=1)                  
plt.title("average pulse")
plt.xlabel("Time past trigger (ms)")
plt.ylabel("Raw counts")
plt.xlim([dt[0], dt[-1]])
plt.ylim(-2000,18000)
plt.grid()        
plt.legend(numpoints=1, frameon=False, loc='best')

ax = plt.subplot(4,2,2)
#np_pulse = np.array(ds.average_pulse)
plt.errorbar(dt[lt:ht], avepulse1[lt:ht], fmt='-', label="Data %s"%dir_pl1, ms=2)
plt.errorbar(dt[lt:ht], avepulse2[lt:ht], fmt='-', label="Data %s"%dir_pl2, ms=2)
plt.errorbar(dt[lt:ht], avepulse3[lt:ht], fmt='-', label="Data %s"%dir_pl3, ms=2)
plt.errorbar(dt[lt:ht], avepulse4[lt:ht], fmt='-', label="Data %s"%dir_pl4, ms=2)                  
plt.title("average pulse")
plt.xlabel("Time past trigger (ms)")
plt.xlim([dt[lt], dt[ht]])
plt.grid()        
plt.legend(numpoints=1, frameon=False, loc='best')

# --- Derivative ---
ac1 = avepulse1 - np.mean(avepulse1)
derivpulse1 =  ac1[1:] - ac1[:-1]
ac2 = avepulse2 - np.mean(avepulse2)
derivpulse2 =  ac2[1:] - ac2[:-1]
ac3 = avepulse3 - np.mean(avepulse3)
derivpulse3 =  ac3[1:] - ac3[:-1]
ac4 = avepulse4 - np.mean(avepulse4)
derivpulse4 =  ac4[1:] - ac4[:-1]
ax = plt.subplot(4,2,3)
plt.errorbar(dt[1:], derivpulse1, fmt='-', label="Data %s"%dir_pl1, ms=1)
plt.errorbar(dt[1:], derivpulse2, fmt='-', label="Data %s"%dir_pl2, ms=1)
plt.errorbar(dt[1:], derivpulse3, fmt='-', label="Data %s"%dir_pl3, ms=1)
plt.errorbar(dt[1:], derivpulse4, fmt='-', label="Data %s"%dir_pl4, ms=1)
plt.title("derivative of average pulse")
plt.xlabel("Time past trigger (ms)")
plt.ylabel("Diff Raw counts")
plt.xlim([dt[0], dt[-1]])
plt.grid()        
plt.legend(numpoints=1, frameon=False, loc='best')

ax = plt.subplot(4,2,4)
plt.errorbar(dt[lt:ht], derivpulse1[lt:ht], fmt='-', label="Data %s"%dir_pl1, ms=2)
plt.errorbar(dt[lt:ht], derivpulse2[lt:ht], fmt='-', label="Data %s"%dir_pl2, ms=2)
plt.errorbar(dt[lt:ht], derivpulse3[lt:ht], fmt='-', label="Data %s"%dir_pl3, ms=2)
plt.errorbar(dt[lt:ht], derivpulse4[lt:ht], fmt='-', label="Data %s"%dir_pl4, ms=2)
plt.title("derivative of average pulse")
plt.xlabel("Time past trigger (ms)")
plt.xlim([dt[lt], dt[ht]])
plt.grid()        
plt.legend(numpoints=1, frameon=False, loc='best')

# --- Template ---
ax = plt.subplot(4,2,5)
plt.errorbar(dt[2:-2], filt1.value, fmt='-', label="Data %s"%dir_pl1, ms=1)
plt.errorbar(dt[2:-2], filt2.value, fmt='-', label="Data %s"%dir_pl2, ms=1)
plt.errorbar(dt[2:-2], filt3.value, fmt='-', label="Data %s"%dir_pl3, ms=1)
plt.errorbar(dt[2:-2], filt4.value, fmt='-', label="Data %s"%dir_pl4, ms=1)
plt.title("template")
plt.xlabel("Time past trigger (ms)")
plt.ylabel("Template")
plt.xlim([dt[0], dt[-1]])
plt.ylim(-0.03,0.03)
plt.grid()        
plt.legend(numpoints=1, frameon=False, loc='best')

ax = plt.subplot(4,2,6)
ac_filt1 = filt1.value - np.mean(filt1.value)
ac_filt1_diff =  ac_filt1[1:] - ac_filt1[:-1]
ac_filt2 = filt2.value - np.mean(filt2.value)
ac_filt2_diff =  ac_filt2[1:] - ac_filt2[:-1]
ac_filt3 = filt3.value - np.mean(filt3.value)
ac_filt3_diff =  ac_filt3[1:] - ac_filt3[:-1]
ac_filt4 = filt4.value - np.mean(filt4.value)
ac_filt4_diff =  ac_filt4[1:] - ac_filt4[:-1]
plt.errorbar(dt[lt:ht], filt1.value[lt:ht], fmt='-', label="Data %s"%dir_pl1, ms=2)
plt.errorbar(dt[lt:ht], filt2.value[lt:ht], fmt='-', label="Data %s"%dir_pl2, ms=2)
plt.errorbar(dt[lt:ht], filt3.value[lt:ht], fmt='-', label="Data %s"%dir_pl3, ms=2)
plt.errorbar(dt[lt:ht], filt4.value[lt:ht], fmt='-', label="Data %s"%dir_pl4, ms=2)
plt.title("template")
plt.xlabel("Time past trigger (ms)")
plt.xlim([dt[lt], dt[ht]])
plt.grid()        
plt.legend(numpoints=1, frameon=False, loc='lower right')


# --- Noise ----
ax = plt.subplot(4,2,7)
freq1 = np.arange(1, 1+len(noisepsd1))*df1
freq2 = np.arange(1, 1+len(noisepsd2))*df2
freq3 = np.arange(1, 1+len(noisepsd3))*df3
freq4 = np.arange(1, 1+len(noisepsd4))*df4
plt.errorbar(freq1, noisepsd1, label="Data %s"%dir_pl1)
plt.errorbar(freq2, noisepsd2, label="Data %s"%dir_pl2)
plt.errorbar(freq3, noisepsd3, label="Data %s"%dir_pl3)
plt.errorbar(freq4, noisepsd4, label="Data %s"%dir_pl4)
ax.set_xlim([freq1[0]*0.9, freq1[-1]*1.1])
units="Counts"
plt.title("noise")
ax.set_ylabel(r'Power Spectral Density ' + '\n' + r'(Counts$^2$/Hz)')
ax.set_xlabel("Frequency (Hz)")
ax.loglog()
plt.grid()        
plt.ylim(5e-4,1e-1)
plt.legend(numpoints=1, frameon=False, loc='best')

ax = plt.subplot(4,2,8)
plt.errorbar(freq1, noisepsd1, label="Data %s"%dir_pl1)
plt.errorbar(freq2, noisepsd2, label="Data %s"%dir_pl2)
plt.errorbar(freq3, noisepsd3, label="Data %s"%dir_pl3)
plt.errorbar(freq4, noisepsd4, label="Data %s"%dir_pl4)
ax.set_xlim([freq1[0]*0.9, freq1[-1]*1.1])
units="Counts"
plt.title("noise")
ax.set_xlabel("Frequency (Hz)")
ax.set_xscale('linear')
ax.set_yscale('linear')        
plt.grid()        
plt.legend(numpoints=1, frameon=False, loc='best')

plt.tight_layout()

# ============================================

plt.figure(figsize=(8.27,11.69))
pmm1 = np.median(premean1[g1])
prm1 = np.median(prerms1[g1])
pmm2 = np.median(premean2[g2])
prm2 = np.median(prerms2[g2])
pmm3 = np.median(premean3[g3])
prm3 = np.median(prerms3[g3])
pmm4 = np.median(premean4[g4])
prm4 = np.median(prerms4[g4])

# ---- pretrig mean hist ----
ax = plt.subplot(2,2,1)
plt.hist(premean1[g1],bins=150,range=(pmm1-50,pmm1+100),normed=True,\
         histtype='step',label="Data %s"%dir_pl1)
plt.hist(premean2[g2],bins=150,range=(pmm1-50,pmm1+100),normed=True,\
         histtype='step',label="Data %s"%dir_pl2)
plt.hist(premean3[g3],bins=150,range=(pmm1-50,pmm1+100),normed=True,\
         histtype='step',label="Data %s"%dir_pl3)
plt.hist(premean4[g4],bins=150,range=(pmm1-50,pmm1+100),normed=True,\
         histtype='step',label="Data %s"%dir_pl4)
ax.set_ylim([0.,0.12])
plt.title("pretrig mean")
ax.set_xlabel("pretrig mean")
ax.set_ylabel("Normalized")
plt.grid()        
plt.legend(numpoints=1, frameon=False, loc='best')

# ---- pretrig mean hist ----
ax = plt.subplot(2,2,2)
plt.hist(premean1[g1]-pmm1,bins=60,range=(-100,200),histtype='step',label="Data %s"%dir_pl1)
plt.hist(premean2[g2]-pmm2,bins=60,range=(-100,200),histtype='step',label="Data %s"%dir_pl2)
plt.hist(premean3[g3]-pmm3,bins=60,range=(-100,200),histtype='step',label="Data %s"%dir_pl3)
plt.hist(premean4[g4]-pmm4,bins=60,range=(-100,200),histtype='step',label="Data %s"%dir_pl4)
plt.title("pretrig mean")
ax.set_xlabel("pretrig mean diff from median")
ax.set_ylabel("Counts / 5 ch")
plt.grid()        
plt.legend(numpoints=1, frameon=False, loc='best')


# ---- pretrig rms hist ----
ax = plt.subplot(2,2,3)
plt.hist(prerms1[g1],bins=80,range=(0,40),normed=True,histtype='step',label="Data %s"%dir_pl1)
plt.hist(prerms2[g2],bins=80,range=(0,40),normed=True,histtype='step',label="Data %s"%dir_pl2)
plt.hist(prerms3[g3],bins=80,range=(0,40),normed=True,histtype='step',label="Data %s"%dir_pl3)
plt.hist(prerms4[g4],bins=80,range=(0,40),normed=True,histtype='step',label="Data %s"%dir_pl4)
plt.title("pretrig rms")
ax.set_ylim([0.,0.6])
ax.set_xlabel("pretrig rms")
ax.set_ylabel("Normalized")
plt.grid()        
plt.legend(numpoints=1, frameon=False, loc='best')

# ---- pretrig rms hist ----
ax = plt.subplot(2,2,4)
plt.hist(prerms1[g1]-prm1,bins=30,range=(-10,20),histtype='step',label="Data %s"%dir_pl1)
plt.hist(prerms2[g2]-prm2,bins=30,range=(-10,20),histtype='step',label="Data %s"%dir_pl2)
plt.hist(prerms3[g3]-prm3,bins=30,range=(-10,20),histtype='step',label="Data %s"%dir_pl3)
plt.hist(prerms4[g4]-prm4,bins=30,range=(-10,20),histtype='step',label="Data %s"%dir_pl4)
plt.title("pretrig rms")
ax.set_xlabel("pretrig rms diff from median")
ax.set_ylabel("Counts / 1 ch")
plt.grid()        
plt.legend(numpoints=1, frameon=False, loc='best')


plt.tight_layout()

plt.figure(figsize=(8.27,11.69))
ax = plt.subplot(2,2,1)
indl = np.logical_and(g3, premean3<np.mean(premean3))
indh = np.logical_and(g3, premean3>np.mean(premean3))
med = np.median(filtvdc3[g3])
plt.hist(filtvdc3[indh],bins=100,range=(med-50,med+50),normed=True,\
         histtype='step',label="high pretrig mean Data %s"%dir_pl3)
plt.hist(filtvdc3[indl],bins=100,range=(med-50,med+50),normed=True,\
         histtype='step',label="low pretrig mean Data %s"%dir_pl3)
ax.set_xlabel("filt_value_dc")
ax.set_ylabel("Norm / 2 ch")
plt.grid()        
plt.legend(numpoints=1, frameon=False, loc='best')

ax = plt.subplot(2,2,2)
indl = np.logical_and(g3, prerms3<np.median(prerms3)+1)
indh = np.logical_and(g3, prerms3>np.median(prerms3)+1)
plt.hist(filtvdc3[indh],bins=50,range=(med-50,med+50),normed=True,\
         histtype='step',label="high pretrig rms Data %s"%dir_pl3)
plt.hist(filtvdc3[indl],bins=50,range=(med-50,med+50),normed=True,\
         histtype='step',label="low pretrig rms Data %s"%dir_pl3)
ax.set_xlabel("filt_value_dc")
ax.set_ylabel("Norm / 2 ch")
plt.grid()        
plt.legend(numpoints=1, frameon=False, loc='best')

ax = plt.subplot(2,2,3)
indl = np.logical_and(g4, premean4<np.mean(premean4))
indh = np.logical_and(g4, premean4>np.mean(premean4))
plt.hist(filtvdc4[indh],bins=50,range=(med-50,med+50),normed=True,\
         histtype='step',label="high pretrig mean Data %s"%dir_pl4)
plt.hist(filtvdc4[indl],bins=50,range=(med-50,med+50),normed=True,\
         histtype='step',label="low pretrig mean Data %s"%dir_pl4)
ax.set_xlabel("filt_value_dc")
ax.set_ylabel("Norm / 2 ch")
plt.grid()        
plt.legend(numpoints=1, frameon=False, loc='best')

ax = plt.subplot(2,2,4)
indl = np.logical_and(g4, prerms4<np.median(prerms4)+1)
indh = np.logical_and(g4, prerms4>np.median(prerms4)+1)
plt.hist(filtvdc4[indh],bins=50,range=(med-50,med+50),normed=True,\
         histtype='step',label="high pretrig rms Data %s"%dir_pl4)
plt.hist(filtvdc4[indl],bins=50,range=(med-50,med+50),normed=True,\
         histtype='step',label="low pretrig rms Data %s"%dir_pl4)
ax.set_xlabel("filt_value_dc")
ax.set_ylabel("Norm / 2 ch")
plt.grid()        
plt.legend(numpoints=1, frameon=False, loc='best')


plt.tight_layout()



# --- save plots ----
dir_save = "/Users/tatsuno/work/heates/data/psites/root"
savedir = dir_save
commands.getoutput('mkdir -p ' + savedir)
print("writing %d plots as pdf to %s"%(len(plt.get_fignums()), savedir))
with PdfPages(path.join(savedir, "Fe55_profiles_chan%03d.pdf"%(ds1.channum))) as pdf:
    for i in plt.get_fignums():
        print("writing plot %d of %d to pdf"%(i, len(plt.get_fignums())))
        pdf.savefig(i)
#print("writing %d plots as png to %s"%(len(plt.get_fignums()), savedir))
#plt.savefig("%s/%s_chan%03d.png"%(savedir,fdate,ds.channum))



pma1 = np.mean(premean1[g1])
pma2 = np.mean(premean2[g2])
pma3 = np.mean(premean3[g3])
pma4 = np.mean(premean4[g4])

pms1 = np.std(premean1[g1])
pms2 = np.std(premean2[g2])
pms3 = np.std(premean3[g3])
pms4 = np.std(premean4[g4])

pra1 = np.mean(prerms1[g1])
pra2 = np.mean(prerms2[g2])
pra3 = np.mean(prerms3[g3])
pra4 = np.mean(prerms4[g4])

prs1 = np.std(prerms1[g1])
prs2 = np.std(prerms2[g2])
prs3 = np.std(prerms3[g3])
prs4 = np.std(prerms4[g4])

print "pretrig mean median  ", pmm1, pmm2, pmm3, pmm4
print "pretrig mean average ", pma1, pma2, pma3, pma4
print "pretrig mean std     ", pms1, pms2, pms3, pms4

print "pretrig rms median   ", prm1, prm2, prm3, prm4
print "pretrig rms average  ", pra1, pra2, pra3, pra4
print "pretrig rms std      ", prs1, prs2, prs3, prs4


h1.close()
h2.close()
h3.close()
h4.close()


# -- end of file
