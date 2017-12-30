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

"""Mainpulates some of the J1713+0747 pulsar data files

Author: Richard M. Prestage
Data:   18 November 2017

"""

import numpy as np
from  GbtRaw import *
import pylab as plt
from scipy.signal import medfilt
from realTimeRfi import *

def main():
    """Show an example of radar interference"""
    chan22 = np.load('chan22.npy')
    nsamples = 100000
    start = 395000
    x_real = chan22[start:start+nsamples,0]
    x_imag = chan22[start:start+nsamples,1]
    P_X = x_real.astype('float')**2 +  x_imag.astype('float')**2

#    mu_LP, SP_loc, WP_loc = impulse_detect3(P_X)

    Fs = 6.25e6
    Ts = 1/Fs * 1e3  # ms

    xpnts = Ts * np.arange(len(P_X)) 
    plt.plot(xpnts,P_X,'b-')
#    plt.plot(xpnts,mu_LP,'r-')
#    plt.plot(xpnts,medfilt(P_X, 101),'g-')
#    plt.plot(Ts * SP_loc, np.zeros_like(SP_loc), 'b.',markersize=10)
#    plt.plot(Ts * WP_loc, np.zeros_like(WP_loc), 'k.',markersize=10)
    plt.xlabel('time in millisec from arbitrary zeropoint')
    plt.ylabel('Power, mean power')
    plt.ylim(-200,10000)
    plt.title('extract from guppi_56465_J1713+0747_0006.0000.raw')
    plt.show()


if __name__ == "__main__":
    main()
    
