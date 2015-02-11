import sys
import time
import datetime
import numpy as np
import pylab as plt
import mass
import cPickle
from os import path
from matplotlib.backends.backend_pdf import PdfPages

import ROOT
from root_numpy import root2array, root2rec, tree2rec

dir_d = "/Users/tatsuno/work/heates/data/psites/root"
fdate = "20141105Q"
ch = 1

timebase=9.6*1e-6
nSamples=1024

# ---- booleans ----
save = False
# ------------------


if __name__=='__main__':
    param = sys.argv
    npar = len(param)
    if (npar==2):
        fdate = str(param[1])
    elif (npar==3):
        fdate = str(param[1])
        ch = str(param[2])
    else:
        print "Error: specify 1.fdate and 2. ch (ex. 20141105P 1)"
        sys.exit(0)


filename = "%s/%s/%s_traces.root"%(dir_d,fdate,fdate)
treename = "at"
#cut = "good==%d && ch==%d"%(1,ch)
cut = "ch==%d"%(ch)

#numpy structured array
arr = root2array(filename,treename=treename,selection=cut)
# good or bad
g = np.where(arr['good']>0)
bad = np.where(arr['good']==0)


# ====================================================================


# plot pulse average vs fv
plt.figure()
plt.plot(arr['pulseave'][g], arr['fv'][g], '.')
plt.plot(arr['pulseave'][bad], arr['fv'][bad], '.')
plt.xlabel("pulse_average")
plt.ylabel("filt_value")
plt.title("chan %d"%ch)

# plot pretrig_mean vs fv
plt.figure()
plt.plot(arr['premean'][g],arr['fv'][g], '.')
plt.plot(arr['premean'][bad],arr['fv'][bad], '.')
plt.xlabel("pretrig mean")
plt.ylabel("filt_value")
plt.title("chan %d"%ch)
plt.grid(which="both",b="on")
plt.minorticks_on()

# plot pretrig_rms vs fv
plt.figure()
plt.plot(arr['prerms'][g],arr['fv'][g], '.')
plt.plot(arr['prerms'][bad],arr['fv'][bad], '.')
plt.xlabel("pretrig rms")
plt.ylabel("filt_value")
plt.title("chan %d"%ch)
plt.grid(which="both",b="on")
plt.minorticks_on()

# plot pretrig_mean vs timestamp
plt.figure()
plt.plot(arr['tstamp'][g],arr['premean'][g], '.')
plt.plot(arr['tstamp'][bad],arr['premean'][bad], '.')
plt.xlabel("timestamp")
plt.ylabel("pretrig mean")
plt.title("chan %d"%ch)
plt.grid(which="both",b="on")
plt.minorticks_on()

# plot pretrig_rms vs timestamp
plt.figure()
plt.plot(arr['tstamp'][g],arr['prerms'][g], '.')
plt.plot(arr['tstamp'][bad],arr['prerms'][bad], '.')
plt.xlabel("timestamp")
plt.ylabel("pretrig rms")
plt.title("chan %d"%ch)
plt.grid(which="both",b="on")
plt.minorticks_on()

# --------------------------------------------------------

# plot first 30 good traces with some condition
plt.figure()
inds = plt.find(np.logical_and(arr['good'], arr['enedc']>6600))[:30]
plt.plot(1e3*timebase*np.arange(nSamples),arr['traces'][inds, :].T)
plt.xlabel("time (ms)")
plt.ylabel("pulse fb signal")
plt.grid(which="both",b="on")
plt.title("chan %d, selected good pulses over 6.6 keV"%ch)
plt.ylim(0,40000)

# plot first 30 good traces with some condition
plt.figure()
inds = plt.find(np.logical_and(arr['good'], np.abs(arr['enedc']-6250)<200))[:30]
plt.plot(1e3*timebase*np.arange(nSamples),arr['traces'][inds, :].T)
plt.xlabel("time (ms)")
plt.ylabel("pulse fb signal")
plt.grid(which="both",b="on")
plt.title("chan %d, selected good pulses 6050<ene<6250"%ch)
plt.ylim(0,40000)

# plot first 30 good traces with some condition
plt.figure()
inds = plt.find(np.logical_and(arr['enedc'][bad]>6050, arr['enedc'][bad]<6250))[:30]
plt.plot(1e3*timebase*np.arange(nSamples),arr['traces'][inds, :].T)
plt.xlabel("time (ms)")
plt.ylabel("pulse fb signal")
plt.grid(which="both",b="on")
plt.title("chan %d, selected bad pulses 6050<ene<6250"%ch)
plt.ylim(0,40000)

# plot first 30 good traces with some condition
plt.figure()
inds = plt.find(np.logical_and(arr['good'], arr['enedc']<4000))[:30]
plt.plot(1e3*timebase*np.arange(nSamples),arr['traces'][inds, :].T)
plt.xlabel("time (ms)")
plt.ylabel("pulse fb signal")
plt.grid(which="both",b="on")
plt.title("chan %d, selected good pulses ene<4keV"%ch)
plt.ylim(0,40000)

# plot first 30 good traces with some condition
plt.figure()
inds = plt.find(np.logical_and(arr['enedc'][bad]>1000, arr['enedc'][bad]<4000))[:30]
plt.plot(1e3*timebase*np.arange(nSamples),arr['traces'][inds, :].T)
plt.xlabel("time (ms)")
plt.ylabel("pulse fb signal")
plt.grid(which="both",b="on")
plt.title("chan %d, selected bad pulses ene<4keV"%ch)
plt.ylim(0,40000)

# plot first 30 good traces with some condition
plt.figure()
inds = plt.find(np.logical_and(arr['good'],np.abs(arr['enedc']-5900)<20.))[:30]
plt.plot(1e3*timebase*np.arange(nSamples),arr['traces'][inds, :].T)
plt.xlabel("time (ms)")
plt.ylabel("pulse fb signal")
plt.grid(which="both",b="on")
plt.title("chan %d, selected good pulses ene=MnKa+-20eV"%ch)
plt.ylim(0,40000)

# plot first 30 good traces with some condition
plt.figure()
inds = plt.find(np.logical_and(arr['good'],np.abs(arr['enedc']-6490)<20.))[:30]
plt.plot(1e3*timebase*np.arange(nSamples),arr['traces'][inds, :].T)
plt.xlabel("time (ms)")
plt.ylabel("pulse fb signal")
plt.grid(which="both",b="on")
plt.title("chan %d, selected good pulses ene=MnKb+-20eV"%ch)
plt.ylim(0,40000)

if save:
    savedir = "%s/%s"%(dir_d,fdate)
    print("writing %d plots as pdf to %s"%(len(plt.get_fignums()), savedir))
    with PdfPages(path.join(savedir, fdate+"_all_figures.pdf")) as pdf:
        for i in plt.get_fignums():
            print("writing plot %d of %d to pdf"%(i, len(plt.get_fignums())))
            pdf.savefig(i)




import commands
params = {'xtick.labelsize': 10, # x ticks
          'ytick.labelsize': 10, # y ticks
          'legend.fontsize': 7
          }
plt.rcParams.update(params)
for j, ds in enumerate(data):
    g = ds.good()

    # get predicted dE
    try:
        if ds.filter is not None:
            rms = ds.filter.variances[filter_name]**0.5
            filterobj = ds.filter
    else:
        rms = ds.hdf5_group['filters/filt_%s'%filter_name].attrs['variance']**0.5
        filterobj = ds.hdf5_group['filters/filt_%s'%filter_name]
        
        print " filterobj = ", filterobj, " values = ", filterobj.value

        rms_fwhm = np.sqrt(np.log(2)*8) # FWHM is this much times the RMS
        v_dv = (1/rms)/rms_fwhm
        print ("Chan %3d filter %-15s Predicted V/dV %6.1f  "
               "Predicted res at %.1f eV: %6.1f eV") % (ds.channum, filter_name, v_dv, std_energy, std_energy/v_dv)
        dE_predict = std_energy/v_dv        

    except Exception, e:
        print "Filter %d can't be used"%j
        print e

    channum = ds.channum
    print "..... channum = ", channum
    if debug: print ">>>>> channum =", channum

    average_pulse = ds.__dict__['average_pulse']
    if debug: print ">>>>> average_pulse =", average_pulse[:]

    # plot average pulse 
        
    noise_psd = ds.__dict__['noise_psd']
    if debug: print ">>>>> noise_psd =", noise_psd[:]

    filename = 'filename'
    filename = ds.__dict__[filename]

    nPulses = 'nPulses'
    nPulses = ds.__dict__[nPulses]

    timestamp_offset = 'timestamp_offset'
    timestamp_offset = ds.__dict__[timestamp_offset]

    print "\n---------------------------------------------------"
    print "len(ds.average_pulse) = ", len(ds.average_pulse)
    print "len(filterobj.value) = ", len(filterobj.value)
    print "ds.nSamples = ", ds.nSamples
    print "ds.nPresamples = ", ds.nPresamples
    print "ds.timebase = ", ds.timebase
    print "-----------------------------------------------------\n"


    dt = (np.arange(ds.nSamples)-ds.nPresamples)*ds.timebase*1e3
        
    startt = 240
    stopt = 310
    cond = np.where ( (dt > dt[startt]) & (dt < dt[stopt]) )

    F = plt.figure(figsize=(10,12))

    ############## average pulse ############################

        ax = plt.subplot(4,2,1)

        titlelabel = "Channel = " + str(channum) +   " nPulses = " + str(nPulses) + " name = " + outputdirname
        plt.figtext(0.3, 0.92, titlelabel)
        plt.figtext(0.7, 0.96, " timestamp_offset = " + str(timestamp_offset), size='small')

#        plt.title(titlelabel) 
        plt.errorbar(dt, ds.average_pulse, fmt='-', label="Chan %d" % ds.channum, ms=1)                    
        plt.xlabel("Time past trigger (ms)")
        plt.ylabel("Raw counts")
        plt.xlim([dt[0], dt[-1]])
        plt.ylim(-2000,18000)
        plt.grid()        

        if use_legend: 
            plt.legend(numpoints=1, frameon=False, loc='best')


        ax = plt.subplot(4,2,2)

        np_pulse = np.array(ds.average_pulse)
        plt.errorbar(dt[cond], np_pulse[cond], fmt='-', label="Chan %d" % ds.channum, ms=2)                    
        plt.xlabel("Time past trigger (ms)")
#        plt.ylabel("Raw counts")
        plt.xlim([dt[startt], dt[stopt]])
#        plt.ylim(-2000,18000)
        plt.grid()        

        if use_legend: 
            plt.legend(numpoints=1, frameon=False, loc='best')

#        plt.errorbar(dt[cond], np_pulse[cond], fmt='ro', label="Chan %d" % ds.channum, ms=2)                    



        ############## diff pulses ############################

        ac = ds.average_pulse - np.mean(ds.average_pulse)
        derivpulse =  ac[1:] - ac[:-1]

        ax = plt.subplot(4,2,3)

        plt.errorbar(dt[1:], derivpulse, fmt='r-', label="Chan %d" % ds.channum, ms=1)                    
 
        plt.xlabel("Time past trigger (ms)")
        plt.ylabel("Diff Raw counts")
        plt.xlim([dt[0], dt[-1]])
#        plt.ylim(-2000,18000)
        plt.grid()        

        if use_legend: 
            plt.legend(numpoints=1, frameon=False, loc='best')


        ax = plt.subplot(4,2,4)

        plt.errorbar(dt[cond], derivpulse[cond], fmt='r-', label="Chan %d" % ds.channum, ms=2)                    
 
        plt.xlabel("Time past trigger (ms)")
#        plt.ylabel("Diff Raw counts")
        plt.xlim([dt[startt], dt[stopt]])
#        plt.ylim(-2000,18000)
        plt.grid()        

        if use_legend: 
            plt.legend(numpoints=1, frameon=False, loc='best')

#        plt.errorbar(dt[cond], derivpulse[cond], fmt='ro', label="Chan %d" % ds.channum, ms=2)                    


        ############## noise PSD ############################

        ax = plt.subplot(4,2,5)

        plt.errorbar(dt[2:-2], filterobj.value, fmt='-', label="Chan %d" % ds.channum, ms=1)                    
        plt.xlabel("Time past trigger (ms)")
        plt.ylabel("Template")
        plt.xlim([dt[0], dt[-1]])
        plt.ylim(-0.03,0.03)
        plt.grid()        

        if use_legend: 
            plt.legend(numpoints=1, frameon=False, loc='best')


        ax = plt.subplot(4,2,6)

        np_filt = np.array(filterobj.value)
        ac_filt = np_filt - np.mean(np_filt)
        ac_filt_diff =  ac_filt[1:] - ac_filt[:-1]

        plt.errorbar(dt[cond], np_filt[cond], fmt='-', label="Chan %d" % ds.channum, ms=2)                    
        plt.errorbar(dt[cond], ac_filt_diff[cond] * 2, fmt='r-', label="DIFF x2 : Chan %d" % ds.channum, ms=2)                    

        plt.xlabel("Time past trigger (ms)")
#        plt.ylabel("Template")
        plt.xlim([dt[startt], dt[stopt]])
#        plt.ylim(-0.03,0.03)
        plt.grid()        

        if use_legend: 
            plt.legend(numpoints=1, frameon=False, loc='best')

#        plt.errorbar(dt[cond], np_filt[cond], fmt='ro', label="Chan %d" % ds.channum, ms=2)                    

        ############## noise PSD ############################

        ax = plt.subplot(4,2,7)

        yvalue = noise_psd 
        df = ds.noise_records.noise_psd.attrs['delta_f']
        freq = np.arange(1, 1+len(yvalue))*df
        print "..... df = ", df 
#        print " freq = ", freq

        plt.errorbar(freq, yvalue, label='TES chan %d' % ds.channum)

        ax.set_xlim([freq[0]*0.9, freq[-1]*1.1])
        units="Counts"
        ax.set_ylabel(r'Power Spectral Density ' + '\n' + r'(Counts$^2$/Hz)')
        ax.set_xlabel("Frequency (Hz)")
        ax.loglog()
        plt.grid()        
        plt.ylim(5e-4,1e-1)

        if use_legend: 
            plt.legend(numpoints=1, frameon=False, loc='best')


        ax = plt.subplot(4,2,8)

        plt.errorbar(freq, yvalue, label='TES chan %d' % ds.channum)

        ax.set_xlim([freq[0]*0.9, freq[-1]*1.1])
        units="Counts"
#        ax.set_ylabel(r'Power Spectral Density ' + '\n' + r'(Counts$^2$/Hz)')
        ax.set_xlabel("Frequency (Hz)")
        ax.set_xscale('linear')
        ax.set_yscale('linear')        
        plt.grid()        
#        plt.ylim(5e-4,1e-1)

        if use_legend: 
            plt.legend(numpoints=1, frameon=False, loc='best')





        ######### save figure #############
        outputfiguredir = outputdirname
        commands.getoutput('mkdir -p ' + outputfiguredir)
        import os
        sfilename = os.path.basename(filename)
        plt.savefig(outputfiguredir + '/' + outputfilename + "_" + str("Chan%03d" % ds.channum) + ".png")
#        plt.savefig(outputfiguredir + '/' + outputdirname + "_" + sfilename.replace('ljh','pdf'))        



