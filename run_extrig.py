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

# --- make exttrig match ROOT file
fname = "%s/%s_trig_match.root"%(dir_p,dir_pl)
a.make_exttrig_match_ROOT(data,fout=fname)



# -- end of file
