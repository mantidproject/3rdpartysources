SET python=%1\python.exe
SET f2py=%1\Scripts\f2py.py

SET files=Back.for chinow.for chosol.for zft.for Deadfit.for dist.for fft.for FileDef.f90 input.for maxent.for Modab.for modamp.for modbak.for move.for Multimaxalpha.for opus.for output.for outspec.for project.for Start.for Tropus.for mantid.for

%python% %f2py% -c %files% --fcompiler=intelvem -m mantid_maxent
