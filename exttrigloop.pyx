#cython: boundscheck=False
#cython: wraparound=False

import time
import datetime

import numpy as np
cimport numpy as np
cimport cython

ctypedef np.int64_t INT_t
ctypedef np.uint64_t UINT_t
ctypedef np.float64_t DOUBLE_t

def external_trig_loop(np.ndarray[DOUBLE_t, ndim=1] tstamp, \
                       np.ndarray[DOUBLE_t, ndim=1] tphase, \
                       np.ndarray[INT_t, ndim=1] stampind, trigtime):
    cdef DOUBLE_t dusec = 50.#micro sec
    cdef DOUBLE_t row_timebase = 0.32*1e-6 # sec
    cdef DOUBLE_t timebase = 9.6*1e-6 # sec
    cdef DOUBLE_t tdiff
    cdef int ir=0
    cdef int tb=0
    cdef int maxtrig=100
    cdef int ip, it, ii
    cdef np.ndarray[UINT_t, ndim=1] tt = np.array(trigtime[ir])
    cdef np.ndarray[DOUBLE_t, ndim=2] tmatch =\
      np.zeros((len(tstamp),maxtrig),dtype=np.float64)
    cdef np.ndarray[DOUBLE_t, ndim=2] diff =\
      np.zeros((len(tstamp),maxtrig),dtype=np.float64)
    cdef DOUBLE_t start = time.time()
    for ip in xrange(len(tstamp)):
        ii=0
        if ip%1000==0:
            print "pulse:%d trig:%d for %.3f sec"%(ip,tb,time.time() - start)
        if len(stampind)!=0 and ip==stampind[ir]:
            print "======= timestamp reset ======="
            ir += 1
            tt = np.array(trigtime[ir])
            print "new trigtime length %d"%(len(tt))
            tb = 0
        if tb==len(tt):continue
        for it in xrange(tb,len(tt),1):
            tdiff = ((tstamp[ip] - tt[it]*row_timebase) + tphase[ip]*timebase)*1e6
            if tdiff<(-1.*dusec):
                tb = it
                break
            elif np.abs(tdiff)<=dusec:
                #print "pulse:%d trig:%d for %.3f sec"%(ip,it,time.time() - start)
                tmatch[ip][ii] = tt[it]*row_timebase
                diff[ip][ii] = tdiff
                ii += 1
                if (ii>=maxtrig):
                    print "Too much number of triggers", maxtrig
                    break
        else:
            tb = len(tt)
            print "End of external trigger loop."
    print "End of pulse loop."
    return tmatch, diff
#    return tmatch



# -- end of file
