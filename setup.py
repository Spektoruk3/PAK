
from setuptools import setup

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

    packages = ['FluidCubeGame'],
    install_requires = ['numpy', 'pygame'],
    )
