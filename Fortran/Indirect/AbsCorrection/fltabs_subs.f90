      FUNCTION F(AMU,T,SEC1,SEC2)
      S=AMU*T*(SEC1-SEC2)
      F=1.0
      IF(S.EQ.0.)then
       F=T
      else
       S=(1-EXP(-S))/S
       F=T*S
      endif
      RETURN
      END
c
      SUBROUTINE open_f(nit,fname)
      character*(*) fname
      logical found
      integer nit
      INQUIRE(FILE=fname,EXIST=found)
      if(found)then
       OPEN(UNIT=nit,FILE=fname,STATUS='OLD',FORM='FORMATTED')
       CLOSE(UNIT=nit,STATUS='DELETE')
      endif
      OPEN(UNIT=nit,FILE=fname,STATUS='NEW',FORM='FORMATTED')
      RETURN
      END
c
