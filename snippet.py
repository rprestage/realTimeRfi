#!/users/rprestag/venv/bin/python

import numpy as np
from  GbtRaw import *
def main():
    """Show an example of radar interference"""
    g = GbtRaw('/home/scratch/rprestag/example.raw')
    data = g.extract(0,3)
    nsamples = 500000

    # Cedric identfied radar intereference starting 24.5M samples
    # from the start of the GUPPI raw data file.
    # correct for the fact that we skipped the first block
    # when creating the little example file.
    # (not sure why we need the factor of 2 with the overlap)
    start = 24500000 - (g.nsamples - 2 * g.overlap)
    end = start + nsamples
    chan22 = (data[22,start:end,:])

    # save to a npy file
    np.save('chan22.npy', chan22)

if __name__ == "__main__":
    main()
    
