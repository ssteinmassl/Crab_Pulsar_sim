#!/usr/bin/env python

import numpy as np
from scipy.stats import expon, norm
from scipy.interpolate import splev, splrep
from matplotlib import pyplot

pyplot.ion()

nsb_rate = 0.4  # GHz p.e.
spe_gain = 10.0  # LSB / p.e. on average

pulse_time, pulse_amplitude = np.loadtxt("flashcam_pulse_shape_fast.dat").T
pulse_dt = np.median(np.diff(pulse_time))
pulse_amplitude *= spe_gain / pulse_amplitude.max()

pulsar_period = 33e-3 / 1e-9  # in nsec
#pulsar_time, pulsar_amplitude = [0, 0.1, 0.3, 0.4, 0.5, 0.9, 1.0], [1, 0, 0, 0.25, 0, 0, 1]  # period, a.u.
pulsar_time, pulsar_amplitude, pulsar_error = np.loadtxt("overallphasogram_profile.dat",unpack = True)
pulsar_mean_err = pulsar_error.mean()
#pulsar_charge_scale = 0.001  # GHz p.e. per a.u.
pulsar_charge_scale = 0.1    # GHz p.e. per a.u.
pulsar_model = splrep(np.array(pulsar_time) * pulsar_period, np.array(pulsar_amplitude) * pulsar_charge_scale, k=1, s=0, per=True)

#t = np.linspace(0, 5 * pulsar_period, num=500)
#pyplot.plot(t, splev(t % pulsar_period, pulsar_model))

time = np.arange(0, 0.001 * pulsar_period,4.0)
trace = np.zeros_like(time)

tmin, tmax, dt = time.min(), time.max(), 4.0
pulse_di = int(np.round(dt / pulse_dt))

t = 0
for delta in expon(scale=1.0 / nsb_rate).rvs(int((tmax - tmin) * nsb_rate + 250)):
    t += delta
    i0 = (t - tmin) / dt
    pulse_i0 = int((i0 % 1) // pulse_dt)
    pulse = pulse_amplitude[pulse_i0::pulse_di]
    pulsar = splev(t % pulsar_period, pulsar_model)
    delta_pulsar = norm(scale = pulsar_mean_err).rvs(1)
    pulse_add =pulse + pulsar + delta_pulsar
    if i0 + pulse.size > trace.size:
        break
    trace[int(i0):int(i0) + pulse.size] += pulse_add

baselines, bl = [], pulse_amplitude.sum() / pulse_di * nsb_rate * dt
for sample in trace:
    if sample > bl:
        bl += 1 / 16
    elif sample < bl:
        bl -= 1/ 16

    baselines.append(bl)

time, trace, baselines = np.array(time), np.array(trace), np.array(baselines)



pyplot.plot(time, trace)
pyplot.plot(time, baselines)
