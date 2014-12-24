import sys
import time
import datetime
import numpy as np
import pylab as plt
import mass
import psidata
import pyroot as a
#import exafs


dir_d = "/Users/tatsuno/work/heates/data"
date_p, ds_p = "20141030","I"
date_n, ds_n = "20141030","F"
dir_pl = "%s_%s"%(date_p,ds_p)
ana = "lowE"

# --- booleans ---
forceNew = False  # update
anado    = False  # do analaysis
newcal   = False  # update calibration
plot     = False  # plot figures for a ch
dump     = True  # dump to ROOT file
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
chans = [i for i in xrange(1,480,2)]
#chans = [1,3]
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

# --- dump to ROOT file
if dump:
    data.set_chan_good(chans)
    fname = "%s/%s_mass.root"%(dir_p,dir_pl)
    a.dump_ROOT(data,fout=fname)

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
