import setuptools
import numpy.distutils.core
import sys

if sys.platform == 'win32':
    extra_link_args = ["-static", "-static-libgfortran", "-static-libgcc"]
else:
    extra_link_args = []

ResNorm = numpy.distutils.core.Extension(name="mantidindirect.bayes.ResNorm",
                                         sources=['mantidindirect/Bayes/ResNorm_main.f90',
                                                  'mantidindirect/Bayes/ResNorm_subs.f90',
                                                  'mantidindirect/Bayes/BlrRes.f90',
                                                  'mantidindirect/Bayes/Bayes.f90',
                                                  'mantidindirect/Bayes/Four.f90',
                                                  'mantidindirect/Bayes/Util.f90'],
                                         extra_link_args=["-static", "-static-libgfortran", "-static-libgcc"])

Quest = numpy.distutils.core.Extension(name="mantidindirect.bayes.Quest",
                                       sources=['mantidindirect/Bayes/Quest_main.f90',
                                                'mantidindirect/Bayes/Quest_subs.f90',
                                                'mantidindirect/Bayes/BlrRes.f90',
                                                'mantidindirect/Bayes/Bayes.f90',
                                                'mantidindirect/Bayes/Four.f90',
                                                'mantidindirect/Bayes/Util.f90',
                                                'mantidindirect/Bayes/Simopt.f90'],
                                       extra_link_args=["-static", "-static-libgfortran", "-static-libgcc"])

QLse = numpy.distutils.core.Extension(name="mantidindirect.bayes.QLse",
                                      sources=['mantidindirect/Bayes/QLse_main.f90',
                                               'mantidindirect/Bayes/QLse_subs.f90',
                                               'mantidindirect/Bayes/BlrRes.f90',
                                               'mantidindirect/Bayes/Bayes.f90',
                                               'mantidindirect/Bayes/Four.f90',
                                               'mantidindirect/Bayes/Util.f90',
                                               'mantidindirect/Bayes/Simopt.f90'],
                                      extra_link_args=["-static", "-static-libgfortran", "-static-libgcc"])

QLres = numpy.distutils.core.Extension(name="mantidindirect.bayes.QLres",
                                       sources=['mantidindirect/Bayes/QLres_main.f90',
                                                'mantidindirect/Bayes/QLres_subs.f90',
                                                'mantidindirect/Bayes/BlrRes.f90',
                                                'mantidindirect/Bayes/Bayes.f90',
                                                'mantidindirect/Bayes/Four.f90',
                                                'mantidindirect/Bayes/Util.f90'],
                                       extra_link_args=["-static", "-static-libgfortran", "-static-libgcc"])

QLdata = numpy.distutils.core.Extension(name="mantidindirect.bayes.QLdata",
                                        sources=['mantidindirect/Bayes/QLdata_main.f90',
                                                 'mantidindirect/Bayes/QLdata_subs.f90',
                                                 'mantidindirect/Bayes/Bayes.f90',
                                                 'mantidindirect/Bayes/Four.f90',
                                                 'mantidindirect/Bayes/Util.f90'],
                                        extra_link_args=["-static", "-static-libgfortran", "-static-libgcc"])

muscat = numpy.distutils.core.Extension(name="mantidindirect.bayes.muscat",
                                        sources=['mantidindirect/AbsCorrection/muscat_data_main.f90',
                                                 'mantidindirect/AbsCorrection/muscat_data.f90',
                                                 'mantidindirect/AbsCorrection/muscat_geom.f90'],
                                        extra_link_args=["-static", "-static-libgfortran", "-static-libgcc"])


from numpy.distutils.command.build_ext import build_ext as _build_ext
class build_ext(_build_ext):
    def finalize_options(self):
        _build_ext.finalize_options(self)
        # If we don't do this on windows, when we do bdist_wheel we wont get a static link
        # this is because it misses the compiler flags to f2py which means it ignores the static flags we try to pass
        if sys.platform == 'win32':
            self.fcompiler = 'gnu95'
            self.compiler = 'mingw32'


numpy.distutils.core.setup(
    name='mantidindirect',
    install_requires=['numpy>=1.17.5'],
    packages=setuptools.find_packages(),
    description='Mantid indirect fortran libraries for Bayes and Absorption correction analysis',
    long_description='These libraries wrap fortran Bayes and Absorption correction libraries using f2py',
    author='Mantid Team',
    ext_modules=[ResNorm, Quest, QLse, QLres, QLdata, muscat],
    author_email="mantid-help@mantidproject.org",
    url='https://www.mantidproject.org',
    version="0.1.5",
    license='BSD',
    cmdclass={'build_ext': build_ext}
)
