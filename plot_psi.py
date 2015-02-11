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
date_p, ds_p = "20141030","I"
date_n, ds_n = "20141030","F"
dir_pl = "%s_%s"%(date_p,ds_p)
ana = "lowE"

# --- booleans ---
forceNew = False  # update
anado    = False  # do analaysis
newcal   = False  # update calibration
resol    = False  # print resolutions
plot     = False  # plot figures for a ch
dump     = False  # dump to ROOT file
# ----------------


if __name__=='__main__':
    param = sys.argv
    npar = len(param)
    if (npar==2):# 20141030_I
        dir_pl = str(param[1])
        date_p = dir_pl.split('_')[0]
        ds_p = dir_pl.split('_')[1]
    elif (npar==3):# 20141030 I
        date_p = str(param[1])
        ds_p = str(param[2])
        dir_pl = "%s_%s"%(date_p,ds_p)
        if param[2]=="True":
            dir_pl = str(param[1])
            date_p = dir_pl.split('_')[0]
            ds_p = dir_pl.split('_')[1]
            forceNew = True
            anado    = True
            newcal   = True
            plot     = False
            dump     = True
    elif (npar==4):# 20141030 I
        date_p = str(param[1])
        ds_p = str(param[2])
        dir_pl = "%s_%s"%(date_p,ds_p)
        if param[3]=="True":
            forceNew = True
            anado    = True
            newcal   = True
            plot     = False
            dump     = True
    else:
        print "Error: specify 1.date and 2.name (ex. 20141030 I)"
        sys.exit(0)

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
chans = [307,325]
exceptch = [83,109,111,115,185,243,285,323,339,349,351,363,397,401,421]
badch = [245,247,275,277]
exceptch.extend(badch)
for i in exceptch:
    if i in chans: chans.remove(i)
if len(chans)==0:
    raise ValueError("zero number of channels selected")
    sys.exit(0)

# --- load data
data = a.load_psi(date_p,ds_p,date_n,ds_n,chs=chans,cut=cuts,ddir=dir_d)

# --- copy this script to mass_output
if "__file__" in locals():
    a.cp_file_to_mass_output(__file__, data.datasets[0].filename)  



# --- analysis
if anado:
    if ana=="Off":
        if dir_pl=="20141031_B":
            other_filter = "%s/%s/ljh_files/%s_%s/%s_%s_mass.hdf5"%(dir_d,date_p,date_p,"C",date_p,"C")
            a.analyze_data_from_other(data,forceNew=forceNew,std_energy=5899,hdf5name=other_filter)
        elif dir_pl=="20141031_D":
            other_filter = "%s/%s/ljh_files/%s_%s/%s_%s_mass.hdf5"%(dir_d,date_p,date_p,"E",date_p,"E")
            a.analyze_data_from_other(data,forceNew=forceNew,std_energy=5899,hdf5name=other_filter)
        elif dir_pl=="20141031_F":
            other_filter = "%s/%s/ljh_files/%s_%s/%s_%s_mass.hdf5"%(dir_d,date_p,date_p,"G",date_p,"G")
            a.analyze_data_from_other(data,forceNew=forceNew,std_energy=5899,hdf5name=other_filter)
        elif dir_pl=="20141031_H":
            other_filter = "%s/%s/ljh_files/%s_%s/%s_%s_mass.hdf5"%(dir_d,date_p,date_p,"G",date_p,"G")
            a.analyze_data_from_other(data,forceNew=forceNew,std_energy=5899,hdf5name=other_filter)
        else: analyze_data(data, forceNew=forceNew, std_energy=5899)
    elif ana=="Mn":
        a.analyze_data_Mn(data,forceNew=forceNew,newcal=newcal)
    elif ana=="Fe":
        a.analyze_data_Fe(data,forceNew=forceNew,newcal=newcal)
    elif ana=="lowE":
        a.analyze_data_lowE(data,forceNew=forceNew,newcal=newcal)
    elif ana=="CrCo":
        a.analyze_data_CrCo(data,forceNew=forceNew,newcal=newcal)
    elif ana=="Se":
        a.analyze_data_Se(data,forceNew=forceNew,newcal=newcal)
    elif ana=="Ge":
        a.analyze_data_Ge(data,forceNew=forceNew,newcal=newcal)
    elif ana=="GaAs":
        a.analyze_data_GaAs(data,forceNew=forceNew,newcal=newcal)
    elif ana=="RbBr":
        a.analyze_data_RbBr(data,forceNew=forceNew,newcal=newcal)
    elif ana=="highE":
        a.analyze_data_highE(data,forceNew=forceNew,newcal=newcal)
    else:
        print "Error: no analysis code"


if resol:
    if ana=="Mn":
        med = a.get_resolutions_Mn(data,nch=len(chans),nmnka=0)
    elif ana=="CrCo":
        #med = a.get_resolutions_Cr(data,nch=len(chans),ncrka=0)
        med = a.get_resolutions_Co(data,nch=len(chans),ncoka=2)
    elif ana=="lowE":
        #med = a.get_resolutions_Mn(data,nch=len(chans),nmnka=1)
        med = a.get_resolutions_Co(data,nch=len(chans),ncoka=5)


        
        
# --- dump to ROOT file
if dump:
    data.set_chan_good(chans)
    fname = "%s/%s_mass.root"%(dir_p,dir_pl)
    a.dump_ROOT(data,fout=fname)




# --------------------------------------------------
params = {'xtick.labelsize': 10, # x ticks
          'ytick.labelsize': 10, # y ticks
          'legend.fontsize': 7
         }
plt.rcParams.update(params)

ds = data.datasets[0];  g = ds.good();
ds1 = data.datasets[1]; g1 = ds1.good();
nsamples=ds.nSamples
presamples=ds.nPresamples
timebase=ds.timebase
dt = (np.arange(nsamples)-presamples)*timebase*1e3
lt=240
ht=310


import h5py
import commands
# --- read HDF5 file ---
h5 = h5py.File("%s/%s_mass.hdf5"%(dir_p,dir_pl))

# ----------------------------------------
ch = h5['chan%d'%ds.channum]
filt = ch['filters']['filt_noconst']
rms = filt.attrs['variance']**0.5
avepulse = ch['average_pulse']
noisepsd = ch['noise_psd']
df = noisepsd.attrs['delta_f']
premean = ch['pretrig_mean']
prerms = ch['pretrig_rms']
timestamp = ch['timestamp']
pulseave = ch['pulse_average']
pulserms = ch['pulse_rms']
filtv = ch['filt_value']
filtvdc = ch['filt_value_dc']
filtvphc = ch['filt_value_phc']
filtvtdc = ch['filt_value_tdc']


ch1 = h5['chan%d'%ds1.channum]
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

# ----------------------------------------

#plt.figure(figsize=(11.69,8.27))
plt.figure(figsize=(8.27,11.69))

# --- Average Pulse ---
ax = plt.subplot(4,2,1)
titlelabel = "Data: %s, target: %s, chan: %d and %d" %(dir_pl,ana,ds.channum,ds1.channum)
plt.figtext(0.04, 0.985, titlelabel, size='small')
#plt.figtext(0.7, 0.96, " timestamp_offset = " + str(timestamp_offset), size='small')
plt.errorbar(dt, avepulse, fmt='-', label="Chan %d"%ds.channum, ms=1)                  
plt.errorbar(dt, avepulse1, fmt='-', label="Chan %d"%ds1.channum, ms=1)                  
plt.title("average pulse")
plt.xlabel("Time past trigger (ms)")
plt.ylabel("Raw counts")
plt.xlim([dt[0], dt[-1]])
plt.ylim(-2000,18000)
plt.grid()        
plt.legend(numpoints=1, frameon=False, loc='best')

ax = plt.subplot(4,2,2)
#np_pulse = np.array(ds.average_pulse)
plt.errorbar(dt[lt:ht], avepulse[lt:ht], fmt='-', label="Chan %d"%ds.channum, ms=2)                    
plt.errorbar(dt[lt:ht], avepulse1[lt:ht], fmt='-', label="Chan %d"%ds1.channum, ms=2)                    
plt.title("average pulse")
plt.xlabel("Time past trigger (ms)")
plt.xlim([dt[lt], dt[ht]])
plt.grid()        
plt.legend(numpoints=1, frameon=False, loc='best')

# --- Derivative ---
ac = avepulse - np.mean(avepulse)
derivpulse =  ac[1:] - ac[:-1]
ac1 = avepulse1 - np.mean(avepulse1)
derivpulse1 =  ac1[1:] - ac1[:-1]
ax = plt.subplot(4,2,3)
plt.errorbar(dt[1:], derivpulse, fmt='-', label="Chan %d"%ds.channum, ms=1)                    
plt.errorbar(dt[1:], derivpulse1, fmt='-', label="Chan %d"%ds1.channum, ms=1)                    
plt.title("derivative of average pulse")
plt.xlabel("Time past trigger (ms)")
plt.ylabel("Diff Raw counts")
plt.xlim([dt[0], dt[-1]])
plt.grid()        
plt.legend(numpoints=1, frameon=False, loc='best')

ax = plt.subplot(4,2,4)
plt.errorbar(dt[lt:ht], derivpulse[lt:ht], fmt='-', label="Chan %d"%ds.channum, ms=2)                    
plt.errorbar(dt[lt:ht], derivpulse1[lt:ht], fmt='-', label="Chan %d"%ds1.channum, ms=2)                    
plt.title("derivative of average pulse")
plt.xlabel("Time past trigger (ms)")
plt.xlim([dt[lt], dt[ht]])
plt.grid()        
plt.legend(numpoints=1, frameon=False, loc='best')

# --- Template ---
ax = plt.subplot(4,2,5)
plt.errorbar(dt[2:-2], filt.value, fmt='-', label="Chan %d"%ds.channum, ms=1)                    
plt.errorbar(dt[2:-2], filt1.value, fmt='-', label="Chan %d"%ds1.channum, ms=1)                    
plt.title("template")
plt.xlabel("Time past trigger (ms)")
plt.ylabel("Template")
plt.xlim([dt[0], dt[-1]])
plt.ylim(-0.03,0.03)
plt.grid()        
plt.legend(numpoints=1, frameon=False, loc='best')

ax = plt.subplot(4,2,6)
ac_filt = filt.value - np.mean(filt.value)
ac_filt_diff =  ac_filt[1:] - ac_filt[:-1]
ac_filt1 = filt1.value - np.mean(filt1.value)
ac_filt1_diff =  ac_filt1[1:] - ac_filt1[:-1]
plt.errorbar(dt[lt:ht], filt.value[lt:ht], fmt='-', label="Chan %d"%ds.channum, ms=2)                    
plt.errorbar(dt[lt:ht], filt1.value[lt:ht], fmt='-', label="Chan %d"%ds1.channum, ms=2)                    
plt.errorbar(dt[lt:ht], ac_filt_diff[lt:ht]*2, fmt='-', label="DIFF x2 : Chan %d"%ds.channum, ms=2)                    
plt.errorbar(dt[lt:ht], ac_filt1_diff[lt:ht]*2, fmt='-', label="DIFF x2 : Chan %d"%ds1.channum, ms=2)                    
plt.title("template")
plt.xlabel("Time past trigger (ms)")
plt.xlim([dt[lt], dt[ht]])
plt.grid()        
plt.legend(numpoints=1, frameon=False, loc='lower right')


# --- Noise ----
ax = plt.subplot(4,2,7)
freq = np.arange(1, 1+len(noisepsd))*df
print "..... df = ", df 
plt.errorbar(freq, noisepsd, label='Chan %d' %ds.channum)
plt.errorbar(freq, noisepsd1, label='Chan %d' %ds1.channum)
ax.set_xlim([freq[0]*0.9, freq[-1]*1.1])
units="Counts"
plt.title("noise")
ax.set_ylabel(r'Power Spectral Density ' + '\n' + r'(Counts$^2$/Hz)')
ax.set_xlabel("Frequency (Hz)")
ax.loglog()
plt.grid()        
plt.ylim(5e-4,1e-1)
plt.legend(numpoints=1, frameon=False, loc='best')

ax = plt.subplot(4,2,8)
plt.errorbar(freq, noisepsd, label='Chan %d' %ds.channum)
plt.errorbar(freq, noisepsd1, label='Chan %d' %ds1.channum)
ax.set_xlim([freq[0]*0.9, freq[-1]*1.1])
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
pmm = np.median(premean[g])
pmm1 = np.median(premean1[g1])
prm = np.median(prerms[g])
prm1 = np.median(prerms1[g1])
plsm = np.median(pulserms[g])
plsm1 = np.median(pulserms1[g1])

# ---- pretrig mean vs time----
ax = plt.subplot(3,2,1)
titlelabel = "Data: %s, target: %s, chan: %d and %d" %(dir_pl,ana,ds.channum,ds1.channum)
plt.figtext(0.04, 0.985, titlelabel, size='small')
plt.errorbar(timestamp[g], premean[g]-pmm, fmt='.', label='Chan %d' %ds.channum)
plt.errorbar(timestamp1[g1], premean1[g1]-pmm1, fmt='.', label='Chan %d' %ds1.channum)
ax.set_xlim([timestamp[0], timestamp[-1]])
plt.title("pretrig_mean vs time")
ax.set_xlabel("time from server start (s)")
ax.set_ylabel("pretrig_mean diff from median")
ax.set_xscale('linear')
ax.set_yscale('linear')        
plt.grid()        
plt.legend(numpoints=1, frameon=False, loc='best')

# ---- pretrig mean hist ----
ax = plt.subplot(3,2,2)
plt.hist(premean[g]-pmm,bins=60,range=(-100,200),histtype='step',label='Chan %d' %ds.channum)
plt.hist(premean1[g1]-pmm1,bins=60,range=(-100,200),histtype='step',label='Chan %d' %ds1.channum)
plt.title("pretrig mean")
ax.set_xlabel("pretrig mean diff from median")
ax.set_ylabel("Counts / 5 ch")
plt.grid()        
plt.legend(numpoints=1, frameon=False, loc='best')


# ---- pretrig rms vs time----
ax = plt.subplot(3,2,3)
plt.errorbar(timestamp[g], prerms[g]-prm, fmt='.', label='Chan %d' %ds.channum)
plt.errorbar(timestamp1[g1], prerms1[g1]-prm1, fmt='.', label='Chan %d' %ds1.channum)
ax.set_xlim([timestamp[0], timestamp[-1]])
plt.title("pretrig rms vs time")
ax.set_xlabel("time from server start (s)")
ax.set_ylabel("pretrig rms diff from median")
ax.set_xscale('linear')
ax.set_yscale('linear')        
plt.grid()        
plt.legend(numpoints=1, frameon=False, loc='best')

# ---- pretrig rms hist ----
ax = plt.subplot(3,2,4)
plt.hist(prerms[g]-prm,bins=30,range=(-10,20),histtype='step',label='Chan %d' %ds.channum)
plt.hist(prerms1[g1]-prm1,bins=30,range=(-10,20),histtype='step',label='Chan %d' %ds1.channum)
plt.title("pretrig rms")
ax.set_xlabel("pretrig rms diff from median")
ax.set_ylabel("Counts / 1 ch")
plt.grid()        
plt.legend(numpoints=1, frameon=False, loc='best')

# ---- pulse rms vs time----
ax = plt.subplot(3,2,5)
plt.errorbar(timestamp[g], pulserms[g]-plsm, fmt='.', label='Chan %d' %ds.channum)
plt.errorbar(timestamp1[g1], pulserms1[g1]-plsm1, fmt='.', label='Chan %d' %ds1.channum)
ax.set_xlim([timestamp[0], timestamp[-1]])
plt.title("pulse rms vs time")
ax.set_xlabel("time from server start (s)")
ax.set_ylabel("pulse rms diff from median")
ax.set_xscale('linear')
ax.set_yscale('linear')        
plt.grid()        
plt.legend(numpoints=1, frameon=False, loc='best')

# ---- pulse rms hist ----
ax = plt.subplot(3,2,6)
plt.hist(pulserms[g]-plsm,bins=60,range=(-150,150),histtype='step',label='Chan %d' %ds.channum)
plt.hist(pulserms1[g1]-plsm1,bins=60,range=(-150,150),histtype='step',label='Chan %d' %ds1.channum)
plt.title("pulse rms")
ax.set_xlabel("pulse rms diff from median")
ax.set_ylabel("Counts / 5 ch")
plt.grid()        
plt.legend(numpoints=1, frameon=False, loc='best')



plt.tight_layout()




# --- save plots ----
dir_save = "/Users/tatsuno/work/heates/data/psites/root"
fdate = "%s%s"%(date_p,ds_p)
savedir = "%s/%s"%(dir_save,fdate)
commands.getoutput('mkdir -p ' + savedir)
print("writing %d plots as pdf to %s"%(len(plt.get_fignums()), savedir))
with PdfPages(path.join(savedir, "%s_chan%03d.pdf"%(fdate,ds.channum))) as pdf:
    for i in plt.get_fignums():
        print("writing plot %d of %d to pdf"%(i, len(plt.get_fignums())))
        pdf.savefig(i)
#print("writing %d plots as png to %s"%(len(plt.get_fignums()), savedir))
#plt.savefig("%s/%s_chan%03d.png"%(savedir,fdate,ds.channum))




h5.close()










# --- ds and plots
'''
sometimes force exit because no calibration data
at ds.calibration['filted value']
'''
if plot:
    ch = chans[0]
    ds = data.channel[ch]
    g = ds.good()
    plt.figure()
    plt.plot(ds.p_pulse_average, ds.p_filt_value, '.')
    plt.plot(ds.p_pulse_average[g], ds.p_filt_value[g], '.')
    plt.xlabel("pulse_average")
    plt.ylabel("p_filt_value")
    plt.title("chan %d"%ch)
    
    plt.figure()
    plt.plot(ds.p_timestamp, '.')
    plt.xlabel("pulse number (arb)")
    plt.ylabel("timestamp (s)")
    plt.title("chan %d"%ch)
    
    plt.figure()
    plt.plot(ds.p_pretrig_mean[g],ds.p_energy[g], '.')
    plt.xlabel("pretrig mean")
    plt.ylabel("energy (eV)")
    plt.title("chan %d"%ch)
    plt.grid(which="both",b="on")
    plt.minorticks_on()

    # plot first 10 good traces
    plt.figure()
    inds = plt.find(g)[:10]
    plt.plot(1e3*ds.timebase*np.arange(ds.nSamples),ds.traces[inds, :].T)
    plt.xlabel("time (ms)")
    plt.ylabel("pulse fb signal")
    plt.grid(which="both",b="on")
    plt.title("chan %d, good pulses"%ch)

    # plot first 10 bad traces
    plt.figure()
    inds = plt.find(~g)[10:20]
    plt.plot(1e3*ds.timebase*np.arange(ds.nSamples),ds.traces[inds, :].T)
    plt.xlabel("time (ms)")
    plt.ylabel("pulse fb signal")
    plt.grid(which="both",b="on")
    plt.title("chan %d, bad pulses"%ch)
    
    # plot first 10 good traces with some condition
    #plt.figure()
    #inds = plt.find(np.logical_and(g, np.abs(ds.p_energy-7000)<50))[:10]
    #plt.plot(1e3*ds.timebase*np.arange(ds.nSamples),ds.traces[inds, :].T)
    #plt.xlabel("time (ms)")
    #plt.ylabel("pulse fb signal")
    #plt.grid(which="both",b="on")
    #plt.title("chan %d, selected good pulses"%ch)

    #calib_crrct = 'p_filt_value_dc'
    #cal = ds.calibration[calib_crrct]
    #fig = mass.calibration.young.diagnose_calibration(cal, hist_plot=True)
    #ds.compare_calibrations()

    data.plot_count_rate()
    #exafs.calibration_summary_compare(data)
    #exafs.calibration_summary(data, calib_crrct)
    ##exafs.pulse_summary(data)
    ##exafs.leftover_phc(data)
    #exafs.cut_vs_time_plot(ds)
    
    # --- save figures on a pdf file
    a.save_plots(data, dir_p, dir_n)



# -- end of file
