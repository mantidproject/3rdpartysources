import setuptools
import numpy.distutils.core


ResNorm = numpy.distutils.core.Extension(name="ResNorm", sources=['mantidindirect/Bayes/ResNorm_main.f90', 'mantidindirect/Bayes/ResNorm_subs.f90',
                                                                  'mantidindirect/Bayes/BlrRes.f90', 'mantidindirect/Bayes/Bayes.f90',
                                                                  'mantidindirect/Bayes/Four.f90', 'mantidindirect/Bayes/Util.f90'])
Quest = numpy.distutils.core.Extension(name="Quest",
                                       sources=['mantidindirect/Bayes/Quest_main.f90',
                                                'mantidindirect/Bayes/Quest_subs.f90',
                                                'mantidindirect/Bayes/BlrRes.f90',
                                                'mantidindirect/Bayes/Bayes.f90',
                                                'mantidindirect/Bayes/Four.f90',
                                                'mantidindirect/Bayes/Util.f90',
                                                'mantidindirect/Bayes/Simopt.f90'])
QLse = numpy.distutils.core.Extension(name="QLse",
                                      sources=['mantidindirect/Bayes/QLse_main.f90',
                                               'mantidindirect/Bayes/QLse_subs.f90',
                                               'mantidindirect/Bayes/BlrRes.f90',
                                               'mantidindirect/Bayes/Bayes.f90',
                                               'mantidindirect/Bayes/Four.f90',
                                               'mantidindirect/Bayes/Util.f90',
                                               'mantidindirect/Bayes/Simopt.f90'])
QLres = numpy.distutils.core.Extension(name="QLres",
                                       sources=['mantidindirect/Bayes/QLres_main.f90',
                                                'mantidindirect/Bayes/QLres_subs.f90',
                                                'mantidindirect/Bayes/BlrRes.f90',
                                                'mantidindirect/Bayes/Bayes.f90',
                                                'mantidindirect/Bayes/Four.f90',
                                                'mantidindirect/Bayes/Util.f90'])

QLdata = numpy.distutils.core.Extension(name="QLdata",
                                        sources=['mantidindirect/Bayes/QLdata_main.f90',
                                                 'mantidindirect/Bayes/QLdata_subs.f90',
                                                 'mantidindirect/Bayes/Bayes.f90',
                                                 'mantidindirect/Bayes/Four.f90',
                                                 'mantidindirect/Bayes/Util.f90'])
muscat = numpy.distutils.core.Extension(name="muscat",
                                        sources=['mantidindirect/AbsCorrection/muscat_data_main.f90',
                                                 'mantidindirect/AbsCorrection/muscat_data.f90',
                                                 'mantidindirect/AbsCorrection/muscat_geom.f90'])
numpy.distutils.core.setup(
    name='mantidindirect',
    ext_modules=[ResNorm, Quest, QLse, QLres, QLdata, muscat],
    description='Mantid indirect fortran libraries for Bayes and Absorption correction analysis',
    author='Mantid Team',
    author_email="mantid-help@mantidproject.org",
    url='https://www.mantidproject.org',
    version="0.1.0",
    license='BSD',
)
