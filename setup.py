"""NAME setup script."""

import os
from setuptools import setup, find_packages

# Get the current version number from inside the module
with open(os.path.join('NAME', 'version.py')) as vf:
    exec(vf.read())

long_description = \
"""
"""

setup(
    name = 'NAME',
    version = __version__,
    description = 'NAME',
    long_description = long_description,
    author = 'Thomas Donoghue',
    author_email = 'tdonoghue.research@gmail.com',
    url = 'https://github.com/voytekresearch/omegamappin',
    packages = find_packages(),
    license = 'Apache License, 2.0',
    classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'License :: OSI Approved :: Apache Software License',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX :: Linux',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7'
    ],
    download_url = 'https://github.com/voytekresearch/omegamappin/releases',
    keywords = ['neuroscience', 'neural oscillations', 'power spectra', '1/f', 'electrophysiology'],
    install_requires = ['numpy', 'scipy', 'pandas', 'seaborn', 'matplotlib', 'fooof'],
    tests_require = ['pytest'],
    extras_require = {}
)
