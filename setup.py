from setuptools import setup, find_packages
import pathlib

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# Read the README file for the long description
long_description = (HERE / "README.md").read_text(encoding="utf-8")

setup(
    name='glorys_plot_tools',
    version='1.0.0',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'cartopy',
        'xarray',
        'cmocean',
        'netCDF4',
        'copernicusmarine',
        'gsw',
    ],
    package_data={
        'glorys_plot_tools': ['*.nc', '*.mat'],
    },
    include_package_data=True,
    author='William Little',
    author_email='willbythesea@gmail.com',
    description='A package for plotting variables and datasets of GLORYS data.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/williamlittle423/glorys_plot_tools',
    license='MIT',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Atmospheric Science',
        'Topic :: Scientific/Engineering :: Visualization',
    ],
    python_requires='>=3.10',
    keywords='GLORYS oceanography plotting visualization',
)
