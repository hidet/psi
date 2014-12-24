import numpy as np
import pylab as plt
import mass
import h5py
from os import path
import shutil
import traceback, sys
import time
import datetime
from scipy.optimize import curve_fit
from sklearn import linear_model
from matplotlib.backends.backend_pdf import PdfPages
import exttrigloop as extl



def cp_file_to_mass_output(fname, ljhfname):
    output_dir = path.dirname(mass.ljh_util.output_basename_from_ljh_fname(ljhfname))
    print(fname, path.join(output_dir, path.split(fname)[-1]))
    shutil.copyfile(fname, path.join(output_dir, path.split(fname)[-1]))
    print "copying %s to %s"%(fname, path.join(output_dir, path.split(fname)[-1]))

    
def save_plots(data, filename, noisename):
    basename = mass.output_basename_from_ljh_fname(data.first_good_dataset.filename)
    dir, fname = path.split(basename)
#    print("writing %d plots as png to %s"%(len(plt.get_fignums()), dir))
    print("writing %d plots as pdf to %s"%(len(plt.get_fignums()), dir))
    with PdfPages(path.join(dir, fname+"_all_figures.pdf")) as pdf:
        for i in plt.get_fignums():
            print("writing plot %d of %d to pdf"%(i, len(plt.get_fignums())))
            pdf.savefig(i)
#            print("writing plot %d of %d as png"%(i, len(plt.get_fignums())))
#            plt.figure(i)
#            plt.savefig(path.join(dir,fname+'figure%d.png') % i, dpi=600)
        d = pdf.infodict()
        d['Title'] = filename
        d['Author'] = noisename
        d['CreationDate'] = datetime.datetime(2009, 11, 13)
        d['ModDate'] = datetime.datetime.today()

        
def load_psi(date_p="20141030", setname_p="A", date_n="20141030", setname_n="B", \
             chs=[1,3], cut="", ddir="/Users/tatsuno/work/heates/data"):
    dir_p = "%s/%s/ljh_files/%s_%s"%(ddir,date_p,date_p,setname_p)
    dir_n = "%s/%s/ljh_files/%s_%s"%(ddir,date_n,date_n,setname_n)
    prefix_p, suffix_p = "%s_%s"%(date_p,setname_p),"ljh"
    prefix_n, suffix_n = "%s_%s"%(date_n,setname_n),"noi"
    pulse_files=["%s/%s_chan%d.%s"%(dir_p, prefix_p, c, suffix_p) for c in chs]
    noise_files=["%s/%s_chan%d.%s"%(dir_n, prefix_n, c, suffix_n) for c in chs]
    data = mass.TESGroup(pulse_files, noise_files)
    data.summarize_data(peak_time_microsec=350.0, forceNew=False)
    data.apply_cuts(cut, forceNew=True)
    print "if you want to check the cuts by gui, run mass.gui.create_cuts(data)"
    return data



def analyze_data(data, forceNew=False, std_energy=5899):
    data.compute_noise_spectra()
    data.avg_pulses_auto_masks() # creates masks and compute average pulses
    data.plot_average_pulses()
    data.compute_filters(f_3db=10000.0,forceNew=forceNew)
    data.summarize_filters(std_energy=std_energy)
    data.filter_data(forceNew=forceNew)
    data.drift_correct(forceNew=forceNew)
    data.phase_correct2014(10, plot=False, forceNew=forceNew, pre_sanitize_p_filt_phase=True)


def analyze_data_from_other(data, forceNew=False, std_energy=5899, hdf5name=""):
    data.compute_noise_spectra()
    data.avg_pulses_auto_masks() # creates masks and compute average pulses
    data.plot_average_pulses()
    data.compute_filters(f_3db=10000.0,forceNew=forceNew)
    data.summarize_filters(std_energy=std_energy)
    #data.filter_data(forceNew=forceNew)
    data.filter_data_from_other(forceNew=forceNew,hdf5name=hdf5name)
    data.drift_correct(forceNew=forceNew)
    data.phase_correct2014(10, plot=False, forceNew=forceNew, pre_sanitize_p_filt_phase=True)


def calib_data(data, forceNew=False, newcal=True, calib=[], excl=[]):
    data.calibrate('p_filt_value_dc', calib, size_related_to_energy_resolution=20.0,
                   excl=excl,forceNew=newcal, max_num_clusters = 18, plot_on_fail=False)
    data.calibrate('p_filt_value_phc', calib, size_related_to_energy_resolution=20.0,
                   excl=excl,forceNew=newcal, max_num_clusters = 18, plot_on_fail=False)
    data.time_drift_correct(forceNew=forceNew)
    data.calibrate('p_filt_value_tdc', calib, size_related_to_energy_resolution=20.0,
                   excl=excl,forceNew=newcal, max_num_clusters = 18, plot_on_fail=False)


def analyze_data_Mn(data, forceNew=False, newcal=True):
    caliblines = ["MnKAlpha", "MnKBeta"]
    excl=[]
    analyze_data(data,forceNew,std_energy=5899)
    calib_data(data,forceNew,newcal,caliblines,excl)

def analyze_data_Fe(data, forceNew=False, newcal=True):
    caliblines = ["FeKAlpha", "FeKBeta"]
    excl=[]
    analyze_data(data,forceNew,std_energy=6404)
    calib_data(data,forceNew,newcal,caliblines,excl)

def analyze_data_Se(data, forceNew=False, newcal=True):
    caliblines = ["SeKAlpha", "SeKBeta"]
    excl=[]
    analyze_data(data,forceNew,std_energy=11222)
#    calib_data(data,forceNew,newcal,caliblines,excl)

def analyze_data_Ge(data, forceNew=False, newcal=True):
    caliblines = ["GeKAlpha", "GeKBeta"]
    excl=[]
    analyze_data(data,forceNew,std_energy=9886)
#    calib_data(data,forceNew,newcal,caliblines,excl)

def analyze_data_GaAs(data, forceNew=False, newcal=True):
    caliblines = ["GaKAlpha", "GaKBeta", "AsKAlpha", "AsKBeta"]
    excl=[]
    analyze_data(data,forceNew,std_energy=9252)
#    calib_data(data,forceNew,newcal,caliblines,excl)

def analyze_data_RbBr(data, forceNew=False, newcal=True):
    caliblines = ["BrKAlpha", "BrKBeta", "RbKAlpha", "RbKBeta"]
    excl=[]
    analyze_data(data,forceNew,std_energy=11924)
#    calib_data(data,forceNew,newcal,caliblines,excl)
        
def analyze_data_CrCo(data, forceNew=False, newcal=True):
    caliblines = ["CrKAlpha", "CrKBeta", "CoKAlpha", "CoKBeta"]
    excl=[]
    analyze_data(data,forceNew,std_energy=5415)
    calib_data(data,forceNew,newcal,caliblines,excl)

def analyze_data_lowE(data, forceNew=False, newcal=True):
    caliblines = ["CrKAlpha", "MnKAlpha", "CrKBeta", "FeKAlpha", "MnKBeta", \
                  "CoKAlpha", "FeKBeta", "CoKBeta", "CuKAlpha", "CuKBeta"]
    excl = ["CrKBeta"]
    analyze_data(data,forceNew,std_energy=5899)
    calib_data(data,forceNew,newcal,caliblines,excl)

def analyze_data_highE(data, forceNew=False, newcal=True):
    caliblines = ["MnKAlpha", "FeKAlpha", "MnKBeta", "FeKBeta", \
                  "GaKAlpha", "GeKAlpha", "GaKBeta", "AsKAplpha", \
                  "GeKBeta", "SeKAlpha", "AsKBeta", "BrKAplpha", \
                  "SeKBeta", "BrKBeta", "RbKAlpha", "RbKBeta"]
    excl = ["MnKBeta", "FeKBeta", "GaKBeta","GeKBeta","SeKBeta","AsKBeta","BrKBeta","RbKBeta"]
    analyze_data(data,forceNew,std_energy=9252)
#    calib_data(data,forceNew,newcal,caliblines,excl)
        



# ------------- dump to ROOT file ---------------

def dump_ROOT(data, fout="tree.root", fopt = "recreate"):
    try:
        import ROOT
    except ImportError:
        raise ValueError('ERROR: cannot import pyROOT')
        return
    ROOT.gROOT.SetBatch(1)
    print "...dump_ROOT start"

    # ----------------------------
    #filter_name = 'noconst'
    #rms_fwhm = np.sqrt(np.log(2)*8) # FWHM is this much times the RMS
    #std_energy=5415. # predicted resolution at std_energy (not importat)
    # ----------------------------
            
    f = ROOT.TFile(fout,fopt)

    cname, ctitle = 'common', 'common tree'
    ct = ROOT.TTree(cname,ctitle)
    bc_ch                 = np.zeros(1,dtype=np.intc)
    bc_col                = np.zeros(1,dtype=np.intc)
    bc_row                = np.zeros(1,dtype=np.intc)     
    bc_ncols              = np.zeros(1,dtype=np.intc) #number_of_columns 8
    bc_nrows              = np.zeros(1,dtype=np.intc) #number_of_rows 30
    bc_npulses            = np.zeros(1,dtype=np.intc) #number_of_pulses
    bc_nsamples           = np.zeros(1,dtype=np.intc) # of sample 1024
    bc_npresamples        = np.zeros(1,dtype=np.intc) # of pre samples 256
    bc_row_timebase       = np.zeros(1,dtype=np.float64)   #row time 0.32us
    bc_timebase           = np.zeros(1,dtype=np.float64)   #frame time 9.6us
    bc_timestamp_offset   = np.zeros(1,dtype=np.float64)   #timestamp offset
    #bc_dE_predict         = np.zeros(1,dtype=np.float64)   
    ct.Branch('ch',              bc_ch               ,'channel number/I')
    ct.Branch('col',             bc_col              ,'column number/I')
    ct.Branch('row',             bc_row              ,'row number/I')
    ct.Branch('ncols',           bc_ncols            ,'n of columns/I')
    ct.Branch('nrows',           bc_nrows            ,'n of rows/I')
    ct.Branch('npulses',         bc_npulses          ,'n of pulses/I')
    ct.Branch('nsample',         bc_nsamples         ,'n of samples/I')
    ct.Branch('npresamples',     bc_npresamples      ,'n of presamples/I')
    ct.Branch('row_timebase',    bc_row_timebase     ,'row timebase/D')
    ct.Branch('timebase',        bc_timebase         ,'frame timebase/D')
    ct.Branch('timestamp_offset',bc_timestamp_offset ,'timestamp offset/D')
    #ct.Branch('dE_predict',      bc_dE_predict       ,'predicted resolution/D')

    start = time.time()
    for ds in data:
        print "channel %d start.... for %.3f"%(ds.channum, (time.time() - start))
        # ---------------------------------------------
        bc_ch[0] = ds.channum
        bc_col[0] = ds.column_number
        bc_row[0] = ds.row_number
        bc_ncols[0] = ds.number_of_columns
        bc_nrows[0] = ds.number_of_rows
        bc_npulses[0] = ds.nPulses
        bc_nsamples[0] = ds.nSamples
        bc_npresamples[0] = ds.nPresamples
        bc_row_timebase[0] = ds.timebase/float(ds.number_of_rows)
        bc_timebase[0] = ds.timebase
        bc_timestamp_offset[0] = ds.timestamp_offset
        #rms = ds.hdf5_group['filters/filt_%s'%filter_name].attrs['variance']**0.5
        #bc_dE_predict[0] = std_energy*rms*rms_fwhm
        ct.Fill()
        # ---------------------------------------------
        
        pname, ptitle = 'chan%d'%ds.channum, 'pulse tree'
        pt = ROOT.TTree(pname,ptitle)
        # --- define branches
        bp_ev = np.zeros(1,dtype=np.intc)
        bp_good = np.zeros(1,dtype=bool)
        bp_filt_phase = np.zeros(1, dtype=np.float64)
        bp_filt_value = np.zeros(1, dtype=np.float64)
        bp_filt_value_dc = np.zeros(1, dtype=np.float64)
        bp_filt_value_phc = np.zeros(1, dtype=np.float64)
        bp_filt_value_tdc = np.zeros(1, dtype=np.float64)
        bp_min_value = np.zeros(1, dtype=np.float64)
        bp_peak_index = np.zeros(1, dtype=np.intc)
        bp_peak_time = np.zeros(1, dtype=np.float64)
        bp_peak_value = np.zeros(1, dtype=np.float64)
        bp_postpeak_deriv = np.zeros(1, dtype=np.float64)
        bp_pretrig_mean = np.zeros(1, dtype=np.float64)
        bp_pretrig_rms = np.zeros(1, dtype=np.float64)
        bp_promptness = np.zeros(1, dtype=np.float64)
        bp_pulse_average = np.zeros(1, dtype=np.float64)
        bp_pulse_rms = np.zeros(1, dtype=np.float64)
        bp_rise_time = np.zeros(1, dtype=np.float64)
        bp_timestamp = np.zeros(1, dtype=np.float64)
        #bp_traces = np.zeros(ds.nSamples, dtype=np.uint16)
        # --- tree branch address
        pt.Branch('ev',            bp_ev, 'event/I')
        pt.Branch('good',          bp_good, 'good/O')
        pt.Branch('filt_phase',    bp_filt_phase, 'filt_phase/D')
        pt.Branch('filt_value',    bp_filt_value, 'filt_value/D')
        pt.Branch('filt_value_dc', bp_filt_value_dc, 'filt_value_dc/D')
        pt.Branch('filt_value_phc',bp_filt_value_phc, 'filt_value_phc/D')
        pt.Branch('filt_value_tdc',bp_filt_value_tdc, 'filt_value_tdc/D')
        pt.Branch('min_value',     bp_min_value, 'min_value/D')
        pt.Branch('peak_index',    bp_peak_index,'peak_index/I')
        pt.Branch('peak_time',     bp_peak_time, 'peak_time/D')
        pt.Branch('peak_value',    bp_peak_value, 'peak_value/D')
        pt.Branch('postpeak_deriv',bp_postpeak_deriv, 'postpeak_deriv/D')
        pt.Branch('pretrig_mean',  bp_pretrig_mean, 'pretrig_mean/D')
        pt.Branch('pretrig_rms',   bp_pretrig_rms, 'pretrig_rms/D')
        pt.Branch('promptness',    bp_promptness, 'promptness/D')
        pt.Branch('pulse_average', bp_pulse_average, 'pulse_average/D')
        pt.Branch('pulse_rms',     bp_pulse_rms, 'pulse_rms/D')
        pt.Branch('rise_time',     bp_rise_time, 'rise_time/D')
        pt.Branch('timestamp',     bp_timestamp, 'timestamp/D')
        #pt.Branch('traces',        bp_traces, 'pulse traces[%d]/s'%ds.nSamples)
        # --- numpy array objects for fast calculation
        ara=np.array(ds.good())
        arb=np.array(ds.p_filt_phase)
        arc=np.array(ds.p_filt_value)
        ard=np.array(ds.p_filt_value_dc)
        are=np.array(ds.p_filt_value_phc)
        arf=np.array(ds.p_filt_value_tdc)
        arg=np.array(ds.p_min_value)
        arh=np.array(ds.p_peak_index)
        ari=np.array(ds.p_peak_time)
        arj=np.array(ds.p_peak_value)
        ark=np.array(ds.p_postpeak_deriv)
        arl=np.array(ds.p_pretrig_mean)
        arm=np.array(ds.p_pretrig_rms)
        arn=np.array(ds.p_promptness)
        aro=np.array(ds.p_pulse_average)
        arp=np.array(ds.p_pulse_rms)
        arq=np.array(ds.p_rise_time)
        arr=np.array(ds.p_timestamp)
        #art=np.array(ds.traces)
        
        for i in xrange(ds.nPulses):
            bp_ev[0]              = i
            bp_good[0]            = ara[i]
            bp_filt_phase[0]      = arb[i]
            bp_filt_value[0]      = arc[i]
            bp_filt_value_dc[0]   = ard[i]
            bp_filt_value_phc[0]  = are[i]
            bp_filt_value_tdc[0]  = arf[i]
            bp_min_value[0]       = arg[i]
            bp_peak_index[0]      = arh[i]
            bp_peak_time[0]       = ari[i]
            bp_peak_value[0]      = arj[i]
            bp_postpeak_deriv[0]  = ark[i]
            bp_pretrig_mean[0]    = arl[i]
            bp_pretrig_rms[0]     = arm[i]
            bp_promptness[0]      = arn[i]
            bp_pulse_average[0]   = aro[i]
            bp_pulse_rms[0]       = arp[i]
            bp_rise_time[0]       = arq[i]
            bp_timestamp[0]       = arr[i]
            #bp_traces             = art[i]

            pt.Fill()
        pt.Write()
    ct.Write()
    f.Close()


# ------- external trigger matching ---------
def get_external_trig_nparray(hdf5name,dsname="trig_times"):
    if not path.exists(hdf5name):
        raise ValueError("Could not find %s"%hdf5name)
        sys.exit(0)
    print "---- reading: %s"%hdf5name
    h5 = h5py.File(hdf5name,'r')
    print "---- making trig_times np.array (np.uint64)"
    trig_times = np.array(h5[dsname], np.uint64)
    print "---- # of trigger: %d"%len(trig_times)
    h5.close()
    return trig_times


def get_reset_index(arr, timebase=0.32*1e-6, nrst=100, tdiff=1e-3):
    '''
    return an array including indecies when server reset happened
    select indices condition (time < 1 msec; tdiff=1e-3)
    this can find triggers very soon after reset 
    '''
    arrst = np.where(arr*timebase<tdiff)
    rstind = np.zeros(nrst,dtype=int)
    if len(arrst[0])==0:# no reset
        rstind = rstind[rstind>0]
        return rstind
    else:
        rstind[0]=arrst[0][0]
        ir=1
        for i, x in enumerate(arrst[0]):
            if i<len(arrst[0])-1 and arrst[0][i+1]>arrst[0][i]+1:
                rstind[ir]=arrst[0][i+1]
                ir = ir+1
        rstind = rstind[rstind>0]
        return rstind


def get_trigtime_list(trig_times, row_timebase=0.32*1e-6,  maxnrst=100, tdiff=1e-3):
    '''
    separate trig_times with resets and return a list of arrays
    '''
    rstind = get_reset_index(trig_times,timebase=row_timebase,nrst=maxnrst,tdiff=tdiff)
    for i,x in enumerate(rstind):
        print "---- external trig reset", i, x, trig_times[x]*row_timebase
    trigtime = []
    if len(rstind)==0:# no reset
        trigtime.append(trig_times[:])
    elif len(rstind)==1:
        trigtime.append(trig_times[:rstind[0]])
        trigtime.append(trig_times[rstind[0]:])
    elif len(rstind)>1:
        trigtime.append(trig_times[:rstind[0]])
        for i in xrange(1,len(rstind),1):
            trigtime.append(trig_times[rstind[i-1]:rstind[i]])
        trigtime.append(trig_times[rstind[len(rstind)-1]:])
    return trigtime


def make_exttrig_match_ROOT(data, fout="hoge_match.root", fopt="recreate"):
    try:
        import ROOT
    except ImportError:
        raise ValueError('ERROR: cannot import pyROOT')
        return
    ROOT.gROOT.SetBatch(1)
    print "...make_exttrig_match_ROOT start"

    start = time.time()
    hdf5name = mass.ljh_util.ljh_get_extern_trig_fname(data.datasets[0].filename)
    trig_times = get_external_trig_nparray(hdf5name)
    print "---- %.3f sec to read hdf5 file"%( time.time() - start )
    maxnrst=100
    resetdiff=1e-3
    row_timebase = data.datasets[0].timebase/float(data.datasets[0].number_of_rows)
    trigtime = get_trigtime_list(trig_times,row_timebase=row_timebase,\
                                 maxnrst=maxnrst,tdiff=resetdiff)
    
    f = ROOT.TFile(fout,fopt)
    hall = ROOT.TH1F("hdiff_all","diff_all",2000,-100,100)

    for ds in data:
        print "---- start chan%d"%ds.channum
        h = ROOT.TH1F("hdiff_ch%d"%ds.channum,"diff_ch%d"%ds.channum,2000,-100,100)
        pname, ptitle = 'chan%d'%ds.channum, 'exttrig match tree'
        pt = ROOT.TTree(pname,ptitle)
        bp_nmatch = np.zeros(1,dtype=np.intc)
        bp_tmatch = np.zeros(maxnrst,dtype=np.float64)
        bp_diff   = np.zeros(maxnrst,dtype=np.float64)
        pt.Branch('nmatch',     bp_nmatch,    'nmatch/I')
        pt.Branch('tmatch',     bp_tmatch,    'tmatch[%d]/D'%maxnrst)
        pt.Branch('diff',       bp_diff,      'diff[%d]/D'%maxnrst)
        timestamp   = np.array(ds.p_timestamp, np.float64)
        filt_phase  = np.array(ds.p_filt_phase, np.float64)
        stmprstind = get_reset_index(timestamp,timebase=1.,nrst=maxnrst,tdiff=1.)
        if len(trigtime)-1 != len(stmprstind):
            print "number of resets mismatching, maybe trigger is strange"
            return
        tmatch, diff = extl.external_trig_loop(timestamp,filt_phase,stmprstind,trigtime)
        for i in xrange(len(tmatch)):# nPulses
            bp_nmatch[0] = len(tmatch[i][np.nonzero(tmatch[i])])
            for j in xrange(bp_nmatch[0]):
                bp_tmatch[j] = tmatch[i][j]
                bp_diff[j] = diff[i][j]
                h.Fill(bp_diff[j])
                hall.Fill(bp_diff[j])
            pt.Fill()
        pt.Write()
        h.Write()
    hall.Write()                    
    hall.DrawCopy()
    f.Close()




