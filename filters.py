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

""" Median Filters

Author: Richard M. Prestage
Data:   18 November 2017

"""

import numpy as np
# import pylab as plt
from scipy.signal import medfilt

def med_hist(hist, nvals, bin_edges):
    """find the median of a histogram"""

    # find histogram bar which contains the median value
    sum = 0
    loc = 0
    for i in range(len(hist)):
        sum += hist[i]
        if sum > nvals / 2:
            loc = i
            break

    # and the median, within this value
    bin_width = bin_edges[1] - bin_edges[0]
    cumm_sum = float(sum-hist[i])
    frac_bin = (float(nvals)/2 - cumm_sum)/hist[i]
    hist_med = bin_edges[loc] + frac_bin * bin_width

    return hist_med

def histfilt(a, kernel=1251, bins=256):
    """median filter using sliding histogram method"""

    nsamples = len(a)
    filt = np.zeros(nsamples)

    # return if kernel is too large
    if nsamples < kernel*2:
        print "kernel is too large"
        return filt

    # we need to use the whole data range for the histogram, 
    # or else it might grow outside the inital values
    hist_range = (a.min(), a.max())

    # make histogram of first "kernel" samples
    hist, bin_edges = np.histogram(a[:kernel], bins=bins,range=hist_range)
    bin_width = bin_edges[1] - bin_edges[0]

    # and find the median
    hist_med = med_hist(hist, kernel, bin_edges)

    # for the first kernel/2 locations, just return this value
    filt[0:kernel/2] = hist_med

    # for the bulk of the array, update the histogram, and use
    # the new median
    for i in range(kernel/2,nsamples-kernel/2):
        # remove the sample falling out of the window
        sample = a[i-kernel/2]
        loc = int((sample-bin_edges[0]-1.0)/bin_width)
        hist[loc] -= 1
        # and add in the new sample
        sample = a[i+kernel/2]
        loc = int((sample-bin_edges[0]-1.0)/bin_width)
        hist[loc] += 1
        # and the new median
        hist_med = med_hist(hist, kernel, bin_edges)
        filt[i] = hist_med

    # for the last kernel/2 locations, just return the latest value
    filt[nsamples-kernel/2:] = hist_med

    # all done
    return filt

def double_histfilt(a, ksmall=625,klarge=12501, bins=256):
    nsamples = a.shape[0]

    print "starting histogram small median filtering"
    hist_filt1 = histfilt(a,ksmall)
    print "done"

    print "starting histogram large median filtering"
    hist_filt2 = histfilt(a,klarge)
    print "done"

    # RFI will only ever increase the power, so take the minimum
    # of the two histogram values
    hist_filt3 = np.zeros(nsamples)
    np.minimum(hist_filt1, hist_filt2, out=hist_filt3)
    return hist_filt3
