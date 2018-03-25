from distutils.core import setup

setup(
    name='pyblock',
    version='0.3',
    author='James Spencer',
    url='http://github.com/jsspencer/pyblock',
    packages=('pyblock', 'pyblock.tests'),
    license='Modified BSD license',
    description='Reblocking analysis tools for correlated data',
    long_description=open('README.rst').read(),
    requires=['numpy', 'pandas (>=0.13)', 'matplotlib'],
)
