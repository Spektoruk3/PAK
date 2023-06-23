
from setuptools import setup
from setuptools import find_packages

version = '0.1.1'

long_description = 'FluidCubeGame'

setup(
    name = 'FluidCubeGame',
    version = version,

    author = 'Leshin Vladimir, Spektoruk Ilya, Krukovskiy Vasiliy',
    author_email = 'i.spektoruk@g.nsu.ru',

    description = 'FluidCubeGame',
    long_description = long_description,

    url = 'https://github.com/Spektoruk3/PAK',

    packages = setuptools.find_packages(),
    install_requires = ['numpy', 'pygame'],
    )
