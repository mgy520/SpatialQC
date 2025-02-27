# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

setup(
    name='SpatialQC',
    version='1.0.3.3',
    description=('Spatial transcriptome quality control, generate HTML QC report, obtain filtered clean data'
    ),
    author='Mao Guangyao',
    url='https://github.com/mgy520/SpatialQC',
    packages=find_packages(),
    include_package_data=True,
    package_data={'SpatialQC': ['data/*.db', 'templates/*.html']},
    entry_points={
        'console_scripts': [
            'SpatialQC = SpatialQC.__main__:main',
        ]
    },
    python_requires='>=3.8,<3.9',
    install_requires=[
            'scanpy>=1.9.8',
            'joblib>=1.4.2',
            'plotly>=5.22.0',
            'scrublet>=0.2.3',
            'anndata2ri>=1.3.1'
            ]
)
