#!/users/rprestag/venv/bin/python

from astropy.io import fits
import numpy as np

def spectroFITS(array, file_name):
    """Writes out array as an image in a FITS file"""
    hdu = fits.PrimaryHDU(array)
    hdu.writeto(file_name)

def main():
    hdulist = fits.open('new.fits')
    hdulist.info()


   # Don't give units or labels, as neither of fv or ds9 seem to do them
   # correctly.
    prihdr = hdulist[0].header
    prihdr['CRPIX1'] = 0.0
    prihdr['CRVAL1'] = 0.0
    prihdr['CDELT1'] = 1.0
    prihdr['CRPIX2'] = 0.0
    prihdr['CRVAL2'] = 1420.0
    prihdr['CDELT2'] = 1.0


    # create passband and time series
    dyn_spec = hdulist[0].data
    print dyn_spec.shape
    bandpass    = np.average(dyn_spec, axis=1)
    timeseries  = np.average(dyn_spec, axis=0)

    # and create new image extensions with these
    bphdu = fits.ImageHDU(bandpass,name='BANDPASS')
    tshdu = fits.ImageHDU(timeseries,name='TIMESERIES')

    # uodate these headers.
    bphdu.header['CRPIX1'] = 0.0
    bphdu.header['CRVAL1'] = 1420.0
    bphdu.header['CDELT1'] = 1.0

    tshdu.header['CRPIX1'] = 0.0
    tshdu.header['CRVAL1'] = 20.0
    tshdu.header['CDELT1'] = 2.0

    hdulist.append(bphdu)
    hdulist.append(tshdu)

    hdulist.writeto('temp.fits')
    hdulist.close()

if __name__ == "__main__":
    main()


