#!/usr/bin/env python3

from setuptools import setup, find_packages
#from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Use README.rst as long_description if available
long_description = ""
try:
    with open(path.join(here, "README.rst"), encoding="utf-8") as f:
        long_description = f.read()
except FileNotFoundError:
    long_description = "A pipeline that uses genome-scale metabolic network modeling to find new targets for antimicrobial drug-design."


setup(
    name = "TARDIS",
    version = "1.0.0",
    description = "A pipeline that uses genome-scale metabolic network modeling to find new targets for antimicrobial drug-design.",
    long_description = long_description,
    long_description_content_type = "text/x-rst",
    license = "CC-BY-4.0",
    author = "Camila Rodrigues Chaves, Artur Duque Rossi, Pedro Henrique Monteiro Torres",
    author_email = "camilachaves@biof.ufrj.br",
    url = "https://github.com/milarchaves/TARDIS",
    packages = find_packages(
        include = [
            "TARDIS", 
            "TARDIS.*"
        ]
    ),
    include_package_data = True,
    package_data = {
        "TARDIS": [
            "Contents/*.html", 
            "Contents/*.svg"
        ]
    },
    keywords = [
        "GEM",
        "targets",
        "antimicrobial",
        "drug-design",
        "network model",
        "metabolic network"
    ],
    classifiers = [
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Creative Commons Attribution 4.0 International",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires = ">=3.9",
    py_modules = [
        "TARDIS.__main__",
        "TARDIS.FindTargets",
        "TARDIS.HomologySearch",
        "TARDIS.Initialise"
    ],
    entry_points = {
        "console_scripts": [
          "tardis = TARDIS.__main__:main"
        ]
    }
)
