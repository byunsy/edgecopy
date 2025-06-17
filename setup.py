import setuptools
from src.edgecopy import __version__, __author__, __license__


setuptools.setup(
    name='edgecopy',
    version=__version__,
    author=__author__,
    license=__license__,
    description='Accurate estimation of copy numbers for duplicated genes using whole-exome sequencing',
    url='https://github.com/byunsy/edgecopy',

    package_dir={'': 'src'},
    packages=setuptools.find_packages(where='src'),
    python_requires='>=3.6',
    include_package_data=True,

    entry_points = dict(console_scripts=['edgecopy=edgecopy.run_all:main']),
    )
