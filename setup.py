"""An open source package for metastable electronic states in molecules.

PyOpenCAP is an open-source package aimed at extending the capabilities of electronic structure packages to 
describe metastable electronic states. 

Currently supported electronic structure packages are:

- Q-Chem
- OpenMolcas
- PySCF
- Psi4

Currently available features:

- Box and Voronoi CAPs
- Trajectory analysis tools

All PyOpenCAP wheels distributed on PyPi are MIT licensed.
"""
DOCLINES = (__doc__ or '').split("\n")
    

import os
import re
import sys
import platform
import subprocess

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion
from os import path

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable,
                      '-DVERSION_INFO=' + self.distribution.get_version()]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
        cmake_args += ['-DBUILD_PYOPENCAP=ON']
        cmake_args += ['-DBUILD_OPENCAP=OFF']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                              self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp)


setup(
    name='pyopencap',
    version='1.2.0',
    author='James Gayvert',
    author_email='jrg444@gmail.com',
    description=DOCLINES[0],
    project_urls = {
        'Source code': "https://github.com/gayverjr/opencap",
        'Documentation': "https://gayverjropencap.readthedocs.io/en/latest/"
    },
    long_description="\n".join(DOCLINES[2:]),
    ext_modules=[CMakeExtension('pyopencap.pyopencap_cpp','opencap')],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
    packages=['pyopencap','pyopencap.analysis'],
    install_requires=['numpy','pandas','scipy'],
)
