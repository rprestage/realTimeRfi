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

def main1():
    """Create a small example GUPPI raw data file

    Copy blocks 1-3  of guppi_56465_J1713+0747_0006.0006.raw 
    to example.raw

"""

    inpath = '/lustre/pulsar/users/rprestag/1713+0747_global/raw/'
    rawfile = 'guppi_56465_J1713+0747_0006.0000.raw'

    g = GbtRaw(inpath + rawfile)
    print "file size in bytes:  ", g.in_size
    print "block size in bytes: ", g.blocsize
    print "number of blocks:    ", g.nblocks
    print "header size in bytes:", g.header_len
    print "hd size in bytes:    ", g.hd_len
    print "number of channels:  ", g.obsnchan
    print "number of samples:   ", g.nsamples

    g.copy('/home/scratch/rprestag/example.raw', 1, 3)

def main2():
    """Show an example of radar interference"""
    g = GbtRaw('/home/scratch/rprestag/example.raw')
    data = g.extract(0,3)

    # Cedric idenfied radar intereference starting 24.5M samples
    # from the start of the GUPPI raw data file.
    # correct for the fact that we skipped the first block
    # when creating the little example file.
    # (not sure why we need the factor of 2 with the overlap)
    start = 24500000 - (g.nsamples - 2 * g.overlap)
    end = start + 500000
    plt.plot(data[22,start:end,2])
    plt.xlabel('samples from arbitrary offset')
    plt.ylabel('Y-pol, real component')
    plt.title('extract from guppi_56465_J1713+0747_0006.0000.raw')
    plt.show()


if __name__ == "__main__":
    main2()
    
