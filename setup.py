from setuptools import find_packages
from numpy.distutils.core import setup

setup(
    name='artemis',
    version='0.0.1',
    description='Molecular allostery Toolkit',
    url='https://github.com/nalsur-veallam/ARTEMIS',
    author='Rislan A. Mallaev',
    author_email='mallaev.ra@phustech.edu',
    license='unkown',

    python_requires=">=3.8",

    packages=find_packages(include=['artemis']),
    install_requires=['numpy', 'matplotlib', 'scipy', 'seaborn', 'pandas'],

    entry_points={
    'console_scripts': [
        'artemis = artemis.main:main'
            ]
    },

    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3.8.10',
    ],
)
