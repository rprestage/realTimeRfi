#!/users/rprestag/venv/bin/python

import numpy as np
import pylab as plt
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

def main():
    ksmall = 625   # 0.1 ms of standard sampling
    kmed   = 3125  # 0.5 ms of standard sampling
    klarge = 12501  # 2 ms of standard sampling
    chan22 = np.load('chan22.npy')
    # extract 50,000 samples
    nsamples = 50000
    start = 395000
    x_real = chan22[start:start+nsamples,0]
    x_imag = chan22[start:start+nsamples,1]
    power = x_real.astype('float')**2 +  x_imag.astype('float')**2

#    print "starting scipy median filtering"
#    sci_filt = medfilt(power,ksmall)
#    print "done"

    print "starting histogram small median filtering"
    hist_filt1 = histfilt(power,ksmall)
    print "done"

    print "starting histogram large median filtering"
    hist_filt2 = histfilt(power,klarge)
    print "done"

    # RFI will only ever increase the power, so take the minimum
    # of the two histogram values
    hist_filt3 = np.zeros(nsamples)
    np.minimum(hist_filt1, hist_filt2, out=hist_filt3)

    plt.plot(power,'g-')
#    plt.plot(sci_filt,'r-')
    plt.plot(hist_filt1,'b-')
    plt.plot(hist_filt2,'r-')
    plt.plot(hist_filt3,'y-')
    plt.show()




if __name__ == "__main__":
    main()
    
