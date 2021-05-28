
# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

import subprocess
import os
import epcy

epcy_version = (
    subprocess.run(["git", "describe", "--tags"], stdout=subprocess.PIPE)
    .stdout.decode("utf-8")
    .strip()
)
assert "." in epcy_version

assert os.path.isfile("epcy/version.py")
with open("epcy/VERSION", "w", encoding="utf-8") as fh:
    fh.write(f"{epcy_version}\n")

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='epcy',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version=epcy_version,

    description='Evaluattion of Predictive CapabilitY for ranking biomarker candidates.',
    long_description=long_description,

    # The project's main homepage.
    url='https://github.com/iric-soft/epcy',

    # Author details
    author='IRIC_bioinfo, Eric Audemard',
    author_email='pipy@iric.ca, eric.audemard@umontreal.ca',

    # Choose your license
    license='MIT',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 4 - Beta',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Intended Audience :: Healthcare Industry',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Software Development',
        'Topic :: Software Development :: Build Tools',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: MIT License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3.6',

        # Others
        'Natural Language :: English',
    ],

    # What does your project relate to?
    keywords='Select predictive indicator',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(exclude=['contrib', 'docs', 'tests']),
    package_data={"epcy": ["VERSION"]},
    include_package_data=True,

    # Alternatively, if you want to distribute just a my_module.py, uncomment
    # this:
    #   py_modules=["my_module"],

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    setup_requires=[
        'cython>=0.29.14',
        'numpy>=1.18.1',
        'pandas>=0.25.3',
        'h5py>=2.10.0',
        'scipy>=1.4.1',
        'scikit-learn>=0.22.1',
        'matplotlib>=3.1.2',
        'numexpr>=2.7.0',
        'seaborn>=0.9.0'
    ],
    install_requires=[
        'cython>=0.29.14',
        'numpy>=1.18.1',
        'pandas>=0.25.3',
        'h5py>=2.10.0',
        'scipy>=1.4.1',
        'scikit-learn>=0.22.1',
        'matplotlib>=3.1.2',
        'numexpr>=2.7.0',
        'seaborn>=0.9.0'
    ],
    python_requires='>3.6',

    # List additional groups of dependencies here (e.g. development
    # dependencies). You can install these using the following syntax,
    # for example:
    # $ pip install -e .[dev,test]

    #TODO
    #extras_require={
    #    'test': ['coverage'],
    #},

    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    # package_data={
    #     'catalog': ['*.fa'],
    # },

    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages. See:
    # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files # noqa
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    # data_files=[('catalog', ['data/catalog'])],

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={
        'console_scripts': [
            'epcy=epcy.epcy:main',
        ],
    },
)
