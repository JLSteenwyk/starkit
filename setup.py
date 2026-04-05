from os import path
from setuptools import setup, find_packages

from starkit.version import __version__

here = path.abspath(path.dirname(__file__))

CLASSIFIERS = [
    "Operating System :: OS Independent",
    "Intended Audience :: Science/Research",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Programming Language :: Python :: 3.12",
    "Programming Language :: Python :: 3.13",
    "Topic :: Scientific/Engineering",
]

REQUIRES = [
    "biopython>=1.82",
    "pyhmmer>=0.10.0",
    "numpy>=1.24.0",
    "jinja2>=3.1.0",
    "mappy>=2.24",
]

setup(
    name="starkit",
    description="Starship prediction in fungal genomes.",
    long_description="StarKIT: Starship prediction in fungal genomes.",
    long_description_content_type="text/markdown",
    author="Jacob L. Steenwyk",
    author_email="jlsteenwyk@gmail.com",
    url="https://github.com/jlsteenwyk/starkit",
    packages=find_packages(),
    classifiers=CLASSIFIERS,
    entry_points={"console_scripts": ["starkit = starkit.starkit:main"]},
    version=__version__,
    include_package_data=True,
    install_requires=REQUIRES,
    python_requires=">=3.9",
)

## push new version to pypi
# rm -rf dist
# python3 setup.py sdist bdist_wheel --universal
# twine upload dist/* -r pypi
