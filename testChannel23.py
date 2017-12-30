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

import sys
import numpy as np
from  GbtRaw import *
from realTimeRfi import *

def make_copy():
    """copy 0006.0000.raw"""
    path    = '/lustre/pulsar/users/rprestag/1713+0747_global/raw/'
    in_file = 'guppi_56465_J1713+0747_0006.0000.raw'
    fitsRoot = '0006.0000'
    g = GbtRaw(path+in_file)
    n = g.get_num_blocks()
    print "the file has", n, "blocks"

    outpath =  '/lustre/pulsar/users/rprestag/1713+0747_global/mit1/'
    g.copy(outpath+in_file, 0, 20)

def mitigate(fileNum):
    path    = '/hyrule/data/users/rprestage/mit3/'
    in_file = 'guppi_56465_J1713+0747_0006.000' + fileNum + '.raw'
    g = GbtRaw(path+in_file, update=True)
    bmax = g.get_num_blocks()
    print "the path is:", path
    print "the file is:", in_file
    print "the file has", bmax, "blocks"

    # loop over all blocks, mitigating channel 23 (channel 22 in file)

    for block in range(bmax):
        print "processing block: ", block
        header, data = g.get_block(block)
        for chan in [22]:
            orig = data[chan,:,:]
            mod  = impulse_mitigate(orig)
            data[chan,:,:] = mod
        g.put_block(header, data, block)


def main():
    fileNum = sys.argv[1]
    mitigate(fileNum)

if __name__ == "__main__":
    main()
    
