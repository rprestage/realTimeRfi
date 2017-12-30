#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import numpy as np
import matplotlib.pyplot as plt

from pulse_detect import pulse_detect

do_plot=True


filename = "/media/cedric/data/data/wave_CASPER/J1713+0747_ch23_ced.dat"

#int8   = np.dtype('int8')
#dt_cpx = np.dtype([('Re',int8),
#                   ('Im',int8),
#                   ])
#dt_sample = np.dtype([('X',dt_cpx),
#                      ('Y',dt_cpx),
#                      ])
#N_samples = os.stat(filename).st_size // dt_sample.itemsize
#dat = np.memmap(filename, dtype=dt_sample, mode ='r', shape=(N_samples))


#P_X = dat['X']['Re'].astype('float')**2 + dat['X']['Im'].astype('float')**2

Fs = 6.25e6
Ts = 1/Fs * 1e3  # ms

start = 24500000
#start = 26300000
#start = 27300000

N_ech = 500000


N_elements = 4 # Xr, Xi, Yr, Yi
#N_samples = os.stat(filename).st_size // N_elements
#dat = np.memmap(filename, dtype='int8', mode ='r', shape=(N_samples, N_elements))
#dat = dat[start:start+N_ech,:]

dat = np.load('chan23.npy')
N_samples = dat.shape[0]

X = dat[:,:2]
Y = dat[:,2:]


P_X = X[:,0].astype('float')**2 + X[:,1].astype('float')**2
SP, WP, mu_LP = pulse_detect(P_X, 26./32, 2**17-1, return_detections=True)

global_mean_power = np.median(mu_LP[N_ech//2::N_ech//50])

SP = np.where(SP==1)[0]
WP = np.where(WP==1)[0]


if do_plot:
  plt.plot(Ts * np.arange(len(P_X)), P_X)
  plt.plot(Ts * np.arange(len(P_X)), mu_LP)
  if len(SP)>0:
    plt.plot(Ts*SP, np.zeros_like(SP), '.r')
  if len(WP)>0:
    plt.plot(Ts*WP, np.zeros_like(WP), '.k')
  plt.xlabel('time (ms)')
  plt.ylabel('Power (arb.)')
  plt.xlim((0,Ts*(N_ech-2)))
  plt.ylim((-global_mean_power//2,40*global_mean_power))
  plt.show()
