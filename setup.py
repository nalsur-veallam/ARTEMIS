import subprocess
from setuptools import setup, find_packages
from setuptools.command.install import install

class InstallCppCommand(install):
    def run(self):
        subprocess.check_call(['make'])
        install.run(self)

setup(
    name='artemis',
    version='0.0.1',
    packages=find_packages(),
    cmdclass={
        'install': InstallCppCommand,
    },
)
