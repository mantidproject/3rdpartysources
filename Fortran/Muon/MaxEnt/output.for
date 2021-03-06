C
C*****************************************************************
C
      SUBROUTINE OUTPUT(P)
C      use FileDef
	common/pulseshape/convolR(8192),CONVOLI(8192)
      common/pulses/npulse,DEF
      common/moments/bmindef
      common/points/npts,ngroups,nhists
      COMMON/RUNDATA/RES,FRAMES,FNORM,IRUNNO,HISTS
      common/flags/FITDEAD,FIXPHASE,FITAMP
      common/ans/ans
      common/amps/amp
      common/sense/phi(96),TAUD(96)
      common/fac/factor,facdef,facfake,ratio
	common/MaxPage/n,f
      COMMON/DETECT/A,B,E,c,d
      REAL A(96),B(96),c(96),d(96),E(8192)
      REAL F(68000),amp(96),asym(96),sum,ao(96)
      INTEGER N,IRUNNO
      CHARACTER*20 FNAME
      character*1 fitamp,ans,an
      fperchan=1./(RES*float(npts)*2.)
!      if ( fixphase.ne.'Y') then
!        write(99,*) ' Save fitted phases ? [N]: '
!        read(5,22) AN
!22      FORMAT(A)
!        IF(AN.EQ.'Y'.OR.AN.EQ.'y') THEN
!          OPEN(15,FILE='PHASE.DAT',STATUS='REPLACE')
!          DO 11 J=1,NGROUPS
!          WRITE(15,*) PHI(J)
!11        CONTINUE
!          CLOSE(15)
!          write(99,*) ' Written phases to phase.dat'
!        ENDIF
!      endif
!      IF (FITDEAD.NE.'Y') GOTO 85
!      write(99,*) ' WRITE NEW DEADTIMES TO DISK ? [N]'
!      read(5,'(A1)') ANS
!      IF (ANS.NE.'Y'.AND.ANS.NE.'y') GOTO 85
!      write(99,*) ' WRITING NEW DEADTIMES TO DISK.'
!      OPEN(16,FILE='TAUD.DAT',STATUS='REPLACE')
!      DO 86 I=1,NGROUPS
!      WRITE(16,*) TAUD(I)
!86    CONTINUE
!      CLOSE(16)
C      OPEN(7,FILE='OUTFILE',STATUS='NEW')    
C      WRITE(7,2)(convolR(I),I=1,N)
C      close(7)
C      OPEN(7,FILE='OUTFILE',STATUS='NEW')    
C      WRITE(7,2)(convolI(I),I=1,N)
C      close(7)
85    WRITE(FNAME,10) IRUNNO
10    FORMAT(I7,'.MAX')
      OPEN(7,FILE=FNAME,STATUS='UNKNOWN')    
      ione=1
      WRITE(7,20) IRUNNO,res,npts,ione,n
20    format(1x,i5,1x,f8.6,1x,i5,1x,i4,1x,i4)
      do 30 i=1,n
30    WRITE(7,2) f(i),float(i-1)*fperchan/135.5e-4
2     FORMAT (E13.4,1x,f10.2)
c

      close(7)
c
      RETURN     	
      END
