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

"""Demo some radar interference in GUPPI data

This routine explores some radar pulses present in channel
23 of file guppi_56465_J1713+0747_0006.0000.raw

Author: Richard M. Prestage
Data:   11 December 2017

"""

import numpy as np
from  GbtRaw import *
import pylab as plt
from matplotlib import rcParams
from realTimeRfi import *

def extract_snippet():
    """Extract a small snippet of the data file"""

    inpath = '/lustre/pulsar/users/rprestag/1713+0747_global/raw/'
    rawfile = 'guppi_56465_J1713+0747_0006.0000.raw'
    g = GbtRaw(inpath + rawfile)
    data = g.extract(0,3)

    # Cedric idenfied radar intereference starting 24.5M samples
    # from the start of the GUPPI raw data file.
    start = 24500000
    end = start + 500000

    # channel 23 on web pages is 22 in disk file
    chan23 = (data[22,start:end,:])

    # save to a npy file
    np.save('chan23.npy', chan23)

def plot_snippet():
    """Plot a small snippet of the data file"""
    chan23 = np.load('chan23.npy')
    x_real = chan23[:,0]
    x_imag = chan23[:,1]
    P_X = x_real.astype('float')**2 +  x_imag.astype('float')**2

    # these come from the header of the data file
    Fs = 6.25e6
    Ts = 1/Fs * 1e3  # ms

    xpnts = Ts * np.arange(len(P_X))
    plt.plot(xpnts, P_X)
    plt.ylim(-200,16000)
    plt.xlabel('time [ms]')
    plt.ylabel('Power [arb. units]')
    plt.title('Extract from guppi_56465_J1713+0747_0006.0000.raw')
    plt.savefig('chan23.png')
    plt.show()


def mitigate1(mode='Viou'):
    """Mitigate some radar pulses using Cedric's and my routine"""

    print "Loading data..."
    chan23 = np.load('chan23.npy')
    nsamples = chan23.shape[0]
    x_real = chan23[:,0]
    x_imag = chan23[:,1]
    P_X = x_real.astype('float')**2 +  x_imag.astype('float')**2

    print "Identifying RFI..."
    if mode == 'Viou':
        SP_loc, WP_loc, mu_LP = pulse_detect(P_X, 26./32, 2**17-1,
                                               return_detections=True)
    elif mode == 'Prestage1':
        mu_LP, SP_loc, WP_loc = impulse_detect1(P_X)
    elif mode == 'Prestage2':
        mu_LP, SP_loc, WP_loc = impulse_detect2(P_X)
    elif mode == 'Prestage3':
        mu_LP, SP_loc, WP_loc = impulse_detect3(P_X)
    elif mode == 'Final':
        mu_LP, SP_loc, WP_loc = impulse_detect(P_X)
    else:
        print "Unknown mode"
    print "done..."


    SP_loc = np.where(SP_loc==1)[0]
    WP_loc = np.where(WP_loc==1)[0]
    SP_frac = 100.0 * len(SP_loc) / nsamples
    WP_frac = 100.0 * len(WP_loc) / nsamples
    print "SP_frac", SP_frac, "WP_frac", WP_frac


    # plot the results

    # set some matplotlib defaults for better appearance
    rcParams.update({'figure.autolayout' : True})
    rcParams.update({'axes.formatter.useoffset' : False})

    Fs = 6.25e6
    Ts = 1/Fs * 1e3  # ms
    xpnts = Ts * np.arange(nsamples)
    plt.plot(xpnts,P_X,'b-')
    plt.plot(xpnts,mu_LP,'r-')
    plt.plot(Ts * SP_loc, np.ones_like(SP_loc)*-200.0, 'r.',markersize=5,
             label='Strong')
    plt.plot(Ts * WP_loc, np.ones_like(WP_loc)*-400.0, 'k.',markersize=5,
             label='Weak')
    plt.ylim(-500,16000)
    plt.xlabel('time [ms]')
    plt.ylabel('Power [arb. units]')
    plt.title('Extract from guppi_56465_J1713+0747_0006.0000.raw')
    plt.text(5,15000, mode)
    plt.legend()
    plt.savefig('chan23_' + mode + '.png')
#    plt.show()

    plt.xlim(40.2, 43.2)
    plt.ylim(-500,4000)
    plt.text(40.5,3500, mode)
    plt.savefig('chan23_' + mode + '_zoom.png')

    plt.show()


def mitigate2():
    """Mitigate some radar pulses and replace with noise"""

    mode = 'mit'

    print "Loading data..."
    chan23 = np.load('chan23.npy')
    nsamples = chan23.shape[0]
    x_real = chan23[:,0]
    x_imag = chan23[:,1]
    P_X = x_real.astype('float')**2 +  x_imag.astype('float')**2

    print "Identifying RFI..."
    mu_LP, SP_loc, WP_loc = impulse_detect(P_X)
    print "done..."

    locs = np.union1d(np.where(SP_loc==1)[0],np.where(WP_loc==1)[0])
    frac = 100.0 * len(locs) / nsamples
    print "frac", frac


    print "replacing RFI"
    r_real = np.copy(x_real)
    r_imag = np.copy(x_imag)
    for i in locs:
        rms = np.sqrt(mu_LP[i] / 2.0)
        r_real[i] = np.rint(np.random.normal(0, rms))
        r_imag[i] = np.rint(np.random.normal(0, rms))

    P_X_mit = r_real.astype('float')**2 +  r_imag.astype('float')**2

    # plot the results

    # set some matplotlib defaults for better appearance
    rcParams.update({'figure.autolayout' : True})
    rcParams.update({'axes.formatter.useoffset' : False})

    Fs = 6.25e6
    Ts = 1/Fs * 1e3  # ms
    xpnts = Ts * np.arange(nsamples)
    plt.plot(xpnts,P_X,'b-')
    plt.plot(xpnts,P_X_mit,'g-')
    plt.plot(xpnts,mu_LP,'r-')
    plt.ylim(-500,16000)
    plt.xlabel('time [ms]')
    plt.ylabel('Power [arb. units]')
    plt.title('Extract from guppi_56465_J1713+0747_0006.0000.raw')
    plt.text(5,15000, mode)
#    plt.legend()
    plt.savefig('chan23_' + mode + '.png')
#    plt.show()

    plt.xlim(40.2, 43.2)
    plt.ylim(-500,4000)
    plt.text(40.5,3500, mode)
    plt.savefig('chan23_' + mode + '_zoom.png')

    plt.xlim(41.3, 41.5)
    plt.ylim(-100,8000)
#    plt.text(40.5,3500, mode)
    plt.savefig('chan23_' + mode + '_zoom2.png')
    plt.show()

    plt.plot(xpnts, x_real,'r-')
    plt.plot(xpnts, r_real, 'b-')
    plt.xlim(40.2, 43.2)
    plt.xlabel('time [ms]')
    plt.ylabel('Real Volts [arb. units]')
    plt.title('Extract from guppi_56465_J1713+0747_0006.0000.raw')
    plt.savefig('chan23_volts.png')
    plt.show()


def mitigate3():
    """Mitigate some radar pulses and replace with noise"""

    mode = 'mit'

    print "Loading data..."
    chan23 = np.load('chan23.npy')
    x_real = chan23[:,0]
    mitData = mitigate(chan23)
    r_real = mitData[:,0]
    plt.plot(x_real,'r-')
    plt.plot(r_real, 'b-')
    plt.xlim(250625, 270000)
    plt.xlabel('samples')
    plt.ylabel('Real Volts [arb. units]')
    plt.title('Extract from guppi_56465_J1713+0747_0006.0000.raw')
    plt.savefig('chan23_volts.png')
    plt.show()





def main():
#    mitigate1(mode='Viou')
#    mitigate1(mode='Prestage1')
#    mitigate1(mode='Prestage2')
#    mitigate1(mode='Prestage3')
#    mitigate1(mode='Final')
    mitigate3()
if __name__ == "__main__":
    main()
    
