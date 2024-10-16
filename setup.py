from setuptools import find_packages
from numpy.distutils.core import setup, Extension
#import pybind11

#ext1 = Extension(name = '_test',
                 #sources = ['artemis/cpp/test.cpp'],
                 #include_dirs=[pybind11.get_include()],
                 #language='c++',
                 #extra_compile_args=['-std=c++11','-Wall','-O3'])

setup(
    name='artemis',
    version='0.0.1',
    description='Molecular communication Toolkit',
    url='https://github.com/nalsur-veallam/ARTEMIS',
    author='Ruslan A. Mallaev',
    author_email='mallaev.ra@phystech.edu',
    license='unkown',

    python_requires=">=3.8",

    packages=find_packages(include=['artemis']),
    install_requires=['numpy', 'matplotlib', 'seaborn', 'pandas', 'scipy', 'biopython', 'tqdm'],

    entry_points={
    'console_scripts': [
        'artemis = artemis.main:main'
            ]
    },

    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'Operating System :: Ubuntu',
        'Programming Language :: Python :: 3.8.10',
    ],
)
