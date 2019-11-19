#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot

#Integrator

trace_read = np.load("trace.npy")
baseline_read = np.load("baseline.npy")
time_read = np.load("time.npy")

tr_bas_diff = trace_read-baseline_read
int_window =128
len_data = int_window*int(tr_bas_diff.size /int_window)
tr_int = np.sum(tr_bas_diff[:len_data].reshape(int_window,int(tr_bas_diff.size /int_window)),axis=0)

time_int = time_read[:len_data][::int_window]
#pyplot.plot(time_int,tr_int)

#fold with period

pulsar_period = 33e-3 / 1e-9  # in nsec
time_fold = time_int % pulsar_period

bin_num = 100
counts,bins_1 = np.histogram(time_fold,bins = bin_num)
counts_weighted,bins_2 = np.histogram(time_fold, weights= tr_int,bins=bin_num)
    
pyplot.hist(bins_2[:-1],bins = bins_2, weights = counts_weighted/counts)
    #ax3.set_yscale("log")
pyplot.ylabel("Mean BL")
pyplot.xlabel("Time [ns]")

#pyplot.plot(time_fold,tr_int)

