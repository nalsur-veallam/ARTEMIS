import subprocess

from setuptools import setup
from setuptools.command.build_py import build_py


class BuildCppCommand(build_py):
    def run(self):
        subprocess.check_call(["make"])
        build_py.run(self)


setup(
    cmdclass={
        "build_py": BuildCppCommand,
    },
)
