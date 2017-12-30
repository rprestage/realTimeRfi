#!/users/rprestag/venv/bin/python

from GbtRaw import *
def main():

    path    = '/lustre/pulsar/users/rprestag/1713+0747_global/raw/'
    in_file = 'guppi_56465_J1713+0747_0006.0000.raw'
    g = GbtRaw(path+in_file)
    g.print_header()

if __name__ == "__main__":
    main()


