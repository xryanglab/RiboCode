"""
A setuptools based setup module.
See:

https://packaging.python.org/en/latest/distributing.html

https://github.com/pypa/sampleproject

"""

# Always prefer setuptools over distutils

from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file

with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()


setup(

    name='RiboCode',

    # Versions should comply with PEP440.  For a discussion on single-sourcing

    # the version across setup.py and the project code, see

    # https://packaging.python.org/en/latest/single_source_version.html

    version='1.2.2',


    #TODO fix name
    description='A package for identifying the translated ORFs using ribosome-profiling data',

    #TODO fix README.rst
    long_description=long_description,

    # The project's main homepage.
    #TODO fix url 
    url='https://github.com/xzt41/RiboCode',

    # Author details
    author='Zhengtao Xiao',

    #TODO fix email
    author_email='xzt13@mails.tsinghua.edu.cn',

    # Choose your license
    #TODO determine a proper license
    license='MIT',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[

        # How mature is this project? Common values are

        #   3 - Alpha

        #   4 - Beta

        #   5 - Production/Stable

        'Development Status :: 2 - Pre-Alpha',

        # Indicate who your project is intended for

        'Intended Audience :: Science/Research',

        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Environment :: Console',




        # Specify the Python versions you support here. In particular, ensure

        # that you indicate whether you support Python 2, Python 3 or both.

        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',

    ],



    # What does your project relate to?
    #TODO set proper keywords
    keywords='ribo-seq ribosome-profiling ORF',



    # You can just specify the packages manually here if your project is

    # simple. Or you can use find_packages().

    packages=find_packages(exclude=['contrib', 'docs', 'tests']),



    # Alternatively, if you want to distribute just a my_module.py, uncomment

    # this:

    #   py_modules=["my_module"],



    # List run-time dependencies here.  These will be installed by pip when

    # your project is installed. For an analysis of "install_requires" vs pip's

    # requirements files see:

    # https://packaging.python.org/en/latest/requirements.html

    install_requires=['pysam>0.8.4','matplotlib','numpy','scipy','pyfasta','biopython','h5py'],


    # List additional groups of dependencies here (e.g. development

    # dependencies). You can install these using the following syntax,

    # for example:

    # $ pip install -e .[dev,data]

    extras_require={

        # 'dev': ['check-manifest'],
        # 'dev': [],
        # 'data': ['coverage'],
        # 'data': [],

    },
    package_data = {
        'gtf_files':  ['data/gtf_files_url.txt'],
        'config_file': ['data/config.txt'],
    },

    # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files # noqa

    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'

    #TODO data file?
    #data_files=[('gtf_files', ['data/gtf_files_url.txt']),('config_file',['data/config.txt'])],

    # To provide executable scripts, use entry points in preference to the

    # "scripts" keyword. Entry points provide cross-platform support and allow

    # pip to create the appropriate form of executable for the target platform.

    entry_points={

        'console_scripts': [
            #name determined the name of cmd line direct call
            'RiboCode=RiboCode.RiboCode:main',
            'prepare_transcripts=RiboCode.prepare_transcripts:main',
            'metaplots=RiboCode.metaplots:main',
            'plot_orf_density=RiboCode.plot_orf_density:main'
        ],

    },

)
