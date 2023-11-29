# For info on how to write a setup.py file, check out the link below
# or ask a friendly neighborhood python programmer!
# https://docs.python.org/3.7/distutils/setupscript.html

from setuptools import setup, find_packages
import glob

setup(
    name = "ullyses",
    version = "0.0.1",
    description = "Create ULLYSES data products",
    author = "ULLYSES STScI Team",
    keywords = ["astronomy"],
    classifiers = ['Programming Language :: Python',
                   'Programming Language :: Python :: 3',
                   'Development Status :: 1 - Planning',
                   'Intended Audience :: Science/Research',
                   'Topic :: Scientific/Engineering :: Astronomy',
                   'Topic :: Scientific/Engineering :: Physics',
                   'Topic :: Software Development :: Libraries :: Python Modules'],
    packages = find_packages(),
    package_dir = {"ullyses": "ullyses"},
    install_requires = ["setuptools",
                        "numpy",
                        "cython",
                        "astropy>=5.0.4",
                        "ullyses_utils",
                        "stistools>=1.4.2",
                        "pandas>2.0",
                        "calcos>=3.4.0",
                        "costools",
                        "matplotlib",
                        "pyyaml",
                        "plotly",
                        "pytest",
                        "hasp>=0.9.5"
                        ],
    dependency_links = ["git+https://github.com/spacetelescope/ullyses-utils@main#egg=ullyses_utils-0.1"],
    )
