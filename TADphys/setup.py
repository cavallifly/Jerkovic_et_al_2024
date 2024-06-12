from setuptools import setup
from distutils.core import setup, Extension

def main():

    setup(
        name         = 'TADphys',
        version      = '0.1',
        author       = 'Marco Di Stefano',
        author_email = 'marco.di-stefano@igh.cnrs.fr',
        packages     = ['tadphys', 'tadphys.modelling'],
        platforms = "OS Independent",
        license = "GPLv3",
        description  = 'TADphys is a Python library that allows to perform biophysical 3D modelling of chromatin regions.',
        long_description = (open("README.rst").read()),
        url          = 'https://github.com/MarcoDiS/tadphys',
        download_url = 'https://github.com/MarcoDiS/tadphys/tarball/master',
    )

if __name__ == '__main__':

    exit(main())
