# setup.py

from setuptools import setup, find_packages

setup(
    name='glorys_plot_utils',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'cartopy',
        'xarray',
        'cmocean',
        'netCDF4',
        'copernicusmarine'
    ],
    package_data={'': ['*.nc', '*.mat']},
    author='William Little',
    author_email='willbythesea@gmail.com',
    description='A package for plotting various variables and datasets of GLORYS data.',
    url='https://github.com/yourusername/my_project',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
    ],
)
