#!/users/rprestag/venv/bin/python

from astropy.io import fits
import numpy as np
import pylab as plt
from matplotlib import rcParams

from GbtRaw import *

def spectroFITS(array, tStart, tRes, fStart, fRes, file_name):
    """Writes out array as an image in a FITS file"""

    # create the dynamic spectrum as the primary image
    hdu = fits.PrimaryHDU(array)

    # add the axes information
    hdu.header['CRPIX1'] = 0.0
    hdu.header['CRVAL1'] = tStart
    hdu.header['CDELT1'] = tRes
    hdu.header['CRPIX2'] = 0.0
    hdu.header['CRVAL2'] = fStart
    hdu.header['CDELT2'] = fRes

    # create the bandpass and timeseries
    bandpass    = np.average(array, axis=1)
    timeseries  = np.average(array, axis=0)

    # and create new image extensions with these
    bphdu = fits.ImageHDU(bandpass,name='BANDPASS')
    tshdu = fits.ImageHDU(timeseries,name='TIMESERIES')
    # uodate these headers.
    bphdu.header['CRPIX1'] = 0.0
    bphdu.header['CRVAL1'] = fStart
    bphdu.header['CDELT1'] = fRes
    tshdu.header['CRPIX1'] = 0.0
    tshdu.header['CRVAL1'] = tStart
    tshdu.header['CDELT1'] = tRes


    hdulist = fits.HDUList([hdu, bphdu, tshdu])
    hdulist.writeto(file_name)

def main():

    rcParams.update({'figure.autolayout' : True})
    rcParams.update({'axes.formatter.useoffset' : False})

    pol = 0
    blocks = 20

    # create 512 time bins of 1024 spectra, each with 159 integrations.
    # there will be a few bytes not used.
    nfreq = 1024
    nspec = 1024
    nint  = 159

    path    = '/lustre/pulsar/users/rprestag/1713+0747_global/raw/'
    in_file = 'guppi_56465_J1713+0747_0006.0000.raw'
    fitsRoot = '0006.0000'
    g = GbtRaw(path+in_file)
    n = g.get_num_blocks()
    print "the file has", n, "blocks"

    # get frequency info, etc, from header
    obsfreq = g.header_dict['OBSFREQ']
    obsbw   = g.header_dict['OBSBW']
    obsnchan= g.header_dict['OBSNCHAN']

    # get the data
    tsData = g.extract(0,blocks, overlap=False)

    # required polarization channels                                       
    sp = 2*pol
    ep = sp+2

    # now loop over all channels
    for chan in range(obsnchan):
        print "Processing channel: ", chan

        # if the bandwidth is -ve, the first channel is at high frequency
        fStart  = obsfreq - float(obsbw)/2.0 + (chan *  float(obsbw) / obsnchan)
        fRes = float(obsbw) / (nfreq * obsnchan)
        tsamp   = abs(float(obsnchan)/obsbw) * 1.0e-06  # MHz to Hz
        tStart  = 0
        tRes    = tsamp * nfreq * nint

        # extract the required channel
        chanData = tsData[chan, :, sp:ep]

        # empty list of power spectra
        spec_list = []

        # do the work
        for s in range(nspec):
            winStart = s * (nfreq * nint)
            accum = np.zeros(nfreq)
            for i in range(nint):
                start = winStart + i * nfreq
                end = start + nfreq
                in_arr = np.zeros((nfreq), dtype=np.complex_)
                in_arr.real = chanData[start:end, 0]
                in_arr.imag = chanData[start:end, 1]
                out_arr = np.fft.fftshift(np.fft.fft(in_arr))
                accum += np.abs(out_arr)**2
            spec_list.append(accum/nint)


        # convert back to numpy array and transpose to desired order
        dyn_spec = np.transpose(np.asarray(spec_list))

        # write it to a fits file
        fitsName = fitsRoot + '.c' + str(chan+1) + '.fits'
        spectroFITS(dyn_spec, tStart, tRes, fStart, fRes, fitsName)    

if __name__ == "__main__":
    main()


