from distutils.core import setup
from setuptools import find_packages

setup(
    name='pyhande',
    version='0.1',
    author='HANDE developers',
    packages=find_packages(include=['pyhande', 'pyhande.*']),
    license='Modified BSD license',
    description='Analysis framework for HANDE calculations',
    long_description=open('README.rst').read(),
    install_requires=['numpy', 'scipy', 'pandas', 'pyblock', 'matplotlib', 'statsmodels'],
)
