from setuptools import setup
import os

setup(
    name='texespy',
    version='0.0',
    packages=['texespy'],
    include_package_data=True,
    package_data={
        "texespy": [
            "data/spice/mk/*.tm",
            "data/spice/kernels_esa/**/*",
        ],
    },
    install_requires=['numpy','matplotlib','scipy','spiceypy','astropy'],
)