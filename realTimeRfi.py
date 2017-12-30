#!/users/rprestag/venv/bin/python

# Copyright (C) 2017 Richard M. Prestage

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

""" RFI Detection Algorithms

Utilities to detect RFI in GBT "raw" data. Eventually, this will
include long-time duration narrow-band RFI, but for now it only
handles impulsive, broadband RFI.

The impulsive RFI detector is based on Cedric Viou's algorithm,
decribed at:

   Dumez-Viou, Cedric, Rodolphe Weber, and Philippe Ravier.
   "Multi-Level Pre-Correlation RFI Flagging for Real-Time 
   Implementation on UniBoard". Journal of Astronomical
   Instrumentation 5.04 (2016): 1641019.

impulse_detect is the final version of the algorithm in use, which
has the choice of IIR filter, or median for the power estimate, and
includes the possibility of a dilation factor. The other routines
are earlier, experimental ones retain for posterity.


Author: Richard M. Prestage
Data:   18 November 2017

"""

import numpy as np
# import pylab as plt
from filters import *

def impulse_mitigate(chan):
    """mitigates a single channel's worth of data (X and Y)."""

    # make a copy to return
    mit = np.copy(chan)

    # mitigate X
    nsamples = chan.shape[0]
    x_real = chan[:,0]
    x_imag = chan[:,1]
    P_X = x_real.astype('float')**2 +  x_imag.astype('float')**2
    mu_LP, SP_loc, WP_loc = impulse_detect(P_X)
    locs = np.union1d(np.where(SP_loc==1)[0],np.where(WP_loc==1)[0])
    for i in locs:
        rms = max(np.sqrt(mu_LP[i] / 2.0),1.0)
        mit[i,0] = np.rint(np.random.normal(0, rms))
        mit[i,1] = np.rint(np.random.normal(0, rms))

    # mitigate Y
    nsamples = chan.shape[0]
    y_real = chan[:,2]
    y_imag = chan[:,3]
    P_Y = y_real.astype('float')**2 +  y_imag.astype('float')**2
    mu_LP, SP_loc, WP_loc = impulse_detect(P_Y)
    locs = np.union1d(np.where(SP_loc==1)[0],np.where(WP_loc==1)[0])
    for i in locs:
        rms = np.sqrt(mu_LP[i] / 2.0)
        mit[i,2] = np.rint(np.random.normal(0, rms))
        mit[i,3] = np.rint(np.random.normal(0, rms))

    return mit


def impulse_detect(P, method='IIR', init=2**15, beta=1.0/2**11, 
                   t_weak=30, t_d_weak=25, t_strong=3, t_d_strong=3, 
                   n_dil=2, lambda_bar=4.0, 
                   lambda_d_bar = 26.0/32.0, correction=1.113110001321999):
    """Detects impulsive RFI in an array of power values.

    Impulsive RFI detector. This version works on successive windows 
    (and so the strong and weak windows must be commensurate), and it
    flags *all* samples in a window as bad, if the threshold is tripped 
    for that window.

    Args:

    P:            array of input powers
    method:       method for calculating smooth power: IIR or MEDIAN
    init:         initializor for low pass IIR filter, default: 128**2 + 128**2
    beta:         Forgetting factor for low pass IIR filter
    t_weak:       number of samples in weak detector window
    t_d_weak:     threshold number of samples to trigger weak detection
    t_strong:     number of samples in strong detector window
    t_d_strong:   threshold number of samples to trigger strong detection
    n_dil:        dilation to add, to catch rising and falling edges
    lambda_bar:   RRP and strong detector threshold level
    lambda_d_bar: weak detector threshold level
    correction:   corrects RRP threshold to strong threshold
    """

    # check that the window sizes are commensurate with each other.
    if t_weak % t_strong != 0:
        print "Window sizes must be commensurate"
        return

    # create the empty low-pass IIR filter; initialize the initial elements,
    # and the strong and weak detection arrays
    nsamples = P.size
    mu_LP    = np.zeros(nsamples)  # low-pass mean power level
    mu_LP[0:n_dil+1] = init
    SP_loc   = np.zeros(nsamples) # strong pulse detection array
    WP_loc   = np.zeros(nsamples) # weak pulse detection array

    # initialize the flagss array for the strong and weak power detectors
    flags_strong = np.zeros(t_strong)
    flags_weak   = np.zeros(t_weak)

    # create the median-filter power values if necessary
    if method == 'MEDIAN':
        mu_LP    = double_histfilt(P)

    # iterate across all samples. Start at location n_dil,in case we need to
    # flag right at the start. in the case of n_dil = 0, start at 1 so we
    # can compare to the previous estimate. End at last integer multiple of 
    # the weak window size.

    if n_dil == 0:
        avail = nsamples -1
    else:
        avail = nsamples - (2 * n_dil)
    start = max(1, n_dil)
    nweak = (avail)/t_weak  # number of weak detector windows
    nstrong = t_weak / t_strong  # number of strong windows in one weak window

    # loop over all weak windows
    for win_weak in range(nweak):
        if (win_weak % 10000) == 0:
            print "Window: ", win_weak, "of: ", nweak

        # loop through this weak window
        for i in range(t_weak):
            loc = start + i
            # update the IIR filter if necessary
            if method == 'IIR':
                if (P[loc] > lambda_bar * mu_LP[loc-1]):
                    mu_LP[loc] = mu_LP[loc-1]
                else:
                    mu_LP[loc] = (beta * P[loc]) + ((1-beta) * mu_LP[loc-1]) 
            # and the weak flags
            if (P[loc] > lambda_d_bar * correction * mu_LP[loc-1]):
                flags_weak[i] = 1
            else:
                flags_weak[i] = 0

        # if treshold is triggered, flag *all* of this weak window as RFI
        if (flags_weak.sum() > t_d_weak):
            WP_loc[start-n_dil:start+t_weak+n_dil] = np.ones(t_weak+(2* n_dil))

        # now loop through the appropriate number of strong windows, using
        # the filtered power values calculated above
        for win_strong in range(nstrong):
            start2 = start + (win_strong * t_strong)
            # calculate flags and threshold for this window
            for i in range(t_strong):
                loc = start2 + i
                if (P[loc] > lambda_bar * mu_LP[loc-1]):
                    flags_strong[i] = 1
                else:
                    flags_strong[i] = 0

            # if treshold is triggered, flag *all* of this strong window as RFI
            if (flags_strong.sum() == t_d_strong):
                SP_loc[start2-n_dil:start2+t_strong+n_dil] = \
                    np.ones(t_strong + (2 * n_dil))

        # end of weak window loop, update start location in power buffer
        start += t_weak

    # and return the results
    return mu_LP, SP_loc, WP_loc


# This is Viou's original code, essentially by me (I added some
# diganostic lines, then removed them again.
def pulse_detect(X, Critere, init, return_detections=False):


  N_Ech = X.shape[0];

  ep = 32./2**16;
  #ep = 0.0863;

  Critere_update = 4
  correction     = 1.113110001321999

  N_screen    = 30;
  N_Threshold = 25;

  mu_LP    = np.zeros(N_Ech);
  mu_LP[0] = init;
  LRS3     = np.zeros(3);   # "000" 
  LRS30    = np.zeros(N_screen);  # "000000000000000000000000000000"

  if return_detections:
    WP     = np.zeros(N_Ech);   # Weak pulse
    SP     = np.zeros(N_Ech);   # Strong pulse

  for i in range(1, N_Ech-1):
    # threshold for mean estimate and strong pulse detection
    if X[i]>Critere_update*mu_LP[i-1]:
      LRS3[1:] = LRS3[:-1]
      LRS3[0] = 1
      diff_in = 0;
    else:
      LRS3[1:] = LRS3[:-1]
      LRS3[0] = 0
      diff_in = X[i] - mu_LP[i-1]
    
    # threshold for weak pulse detection
    if X[i]>correction*Critere*mu_LP[i-1]:
      LRS30[1:] = LRS30[:-1]
      LRS30[0] = 1
    else:
      LRS30[1:] = LRS30[:-1]
      LRS30[0] = 0

    
    mu_LP[i] = mu_LP[i-1] + ep * diff_in;
    
    if return_detections:
      if LRS3.sum() == 3:
        SP[i]=1
      if LRS30.sum() > N_Threshold:
        WP[i]=1;

  if return_detections:
    return SP, WP, mu_LP
  else:
    return mu_LP


def impulse_detect1(P, init=2**17-1, beta=1.0/2**11, t_weak=30, t_d_weak=25, 
                   t_strong=3, t_d_strong=3, lambda_bar=4.0, 
                   lambda_d_bar = 26.0/32.0, correction=1.113110001321999):

    """Detects impulsive RFI in an array of power values.

    This version is based on Cedric's Python routine. I think it has two
    concerns. a) when RFI is detected in a window, *all* samples should
    be flagged, not just the last one. b) it uses a sliding window, so
    every sample is tested multiple times.

    Args:

    P:            array of input powers
    init:         initializor for low pass IIR filter, default: 128**2 + 128**2
    beta:         Forgetting factor for low pass IIR filter
    t_weak:       number of samples in weak detector window
    t_d_weak:     threshold number of samples to trigger weak detection
    t_strong:     number of samples in strong detector window
    t_d_strong:   threshold number of samples to trigger strong detection
    lambda_bar:   RRP and strong detector threshold level
    lambda_d_bar: weak detector threshold level
    correction:   corrects RRP threshold to strong threshold
    """

    # create the empty low-pass IIR filter; initialize the first element,
    # and the locations of strong and weak detection locations
    nsamples = P.size
    mu_LP    = np.zeros(nsamples)  # low-pass mean power level
    mu_LP[0] = init
    SP        = np.zeros(nsamples) # strong pulse locations
    WP        = np.zeros(nsamples) # weak pulse locations

    # initialize the flags array for the strong and weak power detectors
    flags_strong = np.zeros(t_strong)
    flags_weak   = np.zeros(t_weak)

    # iterate across all samples. In Viou's paper, the implication
    # is that samples are treated window by window, but in this 
    # implementation there is simply a sliding window.
    # start at location 1, since we need to compare to the previous estimate

    for i in range(1,nsamples):

        # update low-pass filter and strong detector flag
        if (P[i] > lambda_bar * mu_LP[i-1]):
            flags_strong[1:] = flags_strong[:-1]
            flags_strong[0] = 1
            diff_in = 0
        else:
            flags_strong[1:] = flags_strong[:-1]
            flags_strong[0] = 0
            diff_in = P[i] - mu_LP[i-1]

        # update weak detector flag
        if (P[i] > lambda_d_bar * correction * mu_LP[i-1]):
            flags_weak[1:]   = flags_weak[:-1]
            flags_weak[0] = 1

        else:
            flags_weak[1:]   = flags_weak[:-1]
            flags_weak[0] = 0

        # update low-pass filter
        mu_LP[i] = mu_LP[i-1] + beta * diff_in

        # update the strong and weak detection locations
        if (flags_strong.sum() == t_d_strong):
            SP[i] = 1
        if (flags_weak.sum() > t_d_weak):
            WP[i] = 1

    # and return the results
    return mu_LP, SP, WP


def impulse_detect2(P, init=2**15, beta=1.0/2**11, t_weak=30, t_d_weak=25, 
                   t_strong=3, t_d_strong=3, lambda_bar=4.0, 
                   lambda_d_bar = 26.0/32.0, correction=1.113110001321999):
    """Detects impulsive RFI in an array of power values.

    This is my modified version. It works on successive windows (and
    so the strong and weak windows must be commensurate), and it
    flags *all* samples in a window as bad, if the threshold is tripped 
    for that window.

    Args:

    P:            array of input powers
    init:         initializor for low pass IIR filter, default: 128**2 + 128**2
    beta:         Forgetting factor for low pass IIR filter
    t_weak:       number of samples in weak detector window
    t_d_weak:     threshold number of samples to trigger weak detection
    t_strong:     number of samples in strong detector window
    t_d_strong:   threshold number of samples to trigger strong detection
    lambda_bar:   RRP and strong detector threshold level
    lambda_d_bar: weak detector threshold level
    correction:   corrects RRP threshold to strong threshold
    """

    # check that the window sizes are commensurate with each other.
    if t_weak % t_strong != 0:
        print "Window sizes must be commensurate"
        return

    # create the empty low-pass IIR filter; initialize the first element,
    # and the locations of strong and weak detection locations
    nsamples = P.size
    mu_LP    = np.zeros(nsamples)  # low-pass mean power level
    mu_LP[0] = init
    SP_loc   = np.zeros(nsamples) # strong pulse locations
    WP_loc   = np.zeros(nsamples) # weak pulse locations

    # initialize the flagss array for the strong and weak power detectors
    flags_strong = np.zeros(t_strong)
    flags_weak   = np.zeros(t_weak)

    # iterate across all samples. Start at location 1, since we need to 
    # compare to the previous estimate; end at last integer multiple of 
    # the weak window size.
    start = 1
    nweak = (nsamples-1)/t_weak  # number of weak detector windows
    nstrong = t_weak / t_strong  # number of strong windows in one weak window

    # loop over all weak windows
    for win_weak in range(nweak):

        # loop through this  weak window, update the IIR filter and weak flags
        for i in range(t_weak):
            loc = start + i
            if (P[loc] > lambda_bar * mu_LP[loc-1]):
                mu_LP[loc] = mu_LP[loc-1]
            else:
                mu_LP[loc] = (beta * P[loc]) + ((1-beta) * mu_LP[loc-1]) 

            if (P[loc] > lambda_d_bar * correction * mu_LP[loc-1]):
                flags_weak[i] = 1
            else:
                flags_weak[i] = 0

        # if treshold is triggered, flag *all* of this weak window as RFI
        if (flags_weak.sum() > t_d_weak):
            WP_loc[start:start+t_weak] = np.ones(t_weak)

        # now loop through the appropriate number of strong windows, using
        # the filtered power values calculated above
        for win_strong in range(nstrong):
            start2 = start + (win_strong * t_strong)
            # calculate flags and threshold for this window
            for i in range(t_strong):
                loc = start2 + i
                if (P[loc] > lambda_bar * mu_LP[loc-1]):
                    flags_strong[i] = 1
                else:
                    flags_strong[i] = 0

            # if treshold is triggered, flag *all* of this strong window as RFI
            if (flags_strong.sum() == t_d_strong):
                SP_loc[start2:start2 + t_strong] = np.ones(t_strong)

        # end of weak window loop, update start location in power buffer
        start += t_weak

    # and return the results
    return mu_LP, SP_loc, WP_loc

def impulse_detect3(P, kernel=6251, t_weak=30, t_d_weak=25, 
                   t_strong=3, t_d_strong=3, lambda_bar=4.0, 
                   lambda_d_bar = 26.0/32.0, correction=1.113110001321999):
    """Detects impulsive RFI in an array of power values.

    This third version simply median-filters the power values
    (possible with off-line use).

    Args:

    P:            array of input powers
    kernel:       kernel size for median filter (1ms for standard inputs)
    t_weak:       number of samples in weak detector window
    t_d_weak:     threshold number of samples to trigger weak detection
    t_strong:     number of samples in strong detector window
    t_d_strong:   threshold number of samples to trigger strong detection
    lambda_bar:   RRP and strong detector threshold level
    lambda_d_bar: weak detector threshold level
    correction:   corrects RRP threshold to strong threshold
    """

    # check that the window sizes are commensurate with each other.
    if t_weak % t_strong != 0:
        print "Window sizes must be commensurate"
        return

    # create the median-filter power values
    mu_LP    = double_histfilt(P)
    # and the locations of strong and weak detection locations
    nsamples = P.size
    SP_loc   = np.zeros(nsamples) # strong pulse locations
    WP_loc   = np.zeros(nsamples) # weak pulse locations

    # initialize the flagss array for the strong and weak power detectors
    flags_strong = np.zeros(t_strong)
    flags_weak   = np.zeros(t_weak)

    # iterate across all samples. Start at location 1, since we need to 
    # compare to the previous estimate; end at last integer multiple of 
    # the weak window size.
    start = 1
    nweak = (nsamples-1)/t_weak  # number of weak detector windows
    nstrong = t_weak / t_strong  # number of strong windows in one weak window

    # loop over all weak windows
    for win_weak in range(nweak):

        # loop through this  weak window, update the weak flags
        for i in range(t_weak):
            loc = start + i
            if (P[loc] > lambda_d_bar * correction * mu_LP[loc-1]):
                flags_weak[i] = 1
            else:
                flags_weak[i] = 0

        # if treshold is triggered, flag *all* of this weak window as RFI
        if (flags_weak.sum() > t_d_weak):
            WP_loc[start:start+t_weak] = np.ones(t_weak)

        # now loop through the appropriate number of strong windows
        for win_strong in range(nstrong):
            start2 = start + (win_strong * t_strong)
            # calculate flags and threshold for this window
            for i in range(t_strong):
                loc = start2 + i
                if (P[loc] > lambda_bar * mu_LP[loc-1]):
                    flags_strong[i] = 1
                else:
                    flags_strong[i] = 0

            # if treshold is triggered, flag *all* of this strong window as RFI
            if (flags_strong.sum() == t_d_strong):
                SP_loc[start2:start2 + t_strong] = np.ones(t_strong)

        # end of weak window loop, update start location in power buffer
        start += t_weak

    # and return the results
    return mu_LP, SP_loc, WP_loc
