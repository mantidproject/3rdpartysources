MODULE FileDef


LOGICAL(kind=4) PAR_FILE_INP, DATA_FILE_INP, FIT_FILE_INP
CHARACTER(LEN=255) PAR_INP, DATA_INP, cdir, FIT_INP,data_inp_par, wdr
!REAL TAUD(32)

END MODULE FileDef


MODULE MaxDef

character*255 datafilename,readdata
logical onepulse, varydt, fixphasein, retlog
integer nptsin, zeroch, istgoodch, ierr, debug
real deflevel, Sigmaloss,MinFldPlot,MaxFldPlot
character*10 temp, fieldchar
integer runit

end module maxdef

module SaveDef

Real mom1(500),mom2(500),mom3(500),SaveTemp(500),SaveField(500)
integer numsaved
end module SaveDef

