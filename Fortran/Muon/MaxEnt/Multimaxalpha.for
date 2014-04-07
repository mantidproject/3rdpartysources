cc this is multimax.for - jan 96 multiple numor reads and currant bun 
c  found the asymmetry bug 9/1/96: asymmetry with free phases was wrong
c  link with multibunrd
C
C**********************************************************************
C          RAL VERSION OF MAXENT - MAIN PROGRAM : reads runs itself
C**********************************************************************
C      
c      PROGRAM multimax
      subroutine multimax


      use maxdef
     	use dflib
	use dflogm


      REAL DATUM(68000),SIGMA(68000),F(68000),BASE(68000),
     +CONVOLR(8192),CONVOLI(8192),HISTS(64)
      INTEGER P,GROUP(64)
      CHARACTER*1 ANS,AN,FITAMP,firstgo
      COMMON/FILE/NAME
      COMMON/DETECT/A,B,E,c,d
      REAL A(64),B(64),c(64),d(64),E(8192)
      common/sense/phi(64),TAUD(64)
      REAL TAUD
      REAL CORR(68000),DATT(68000)
      logical FITDEAD,FIXPHASE
	character*2046 str
      common/fac/factor,facdef,facfake,ratio
      common/channels/itzero,ipulse1,ipulse2,i1stgood,itotal
      common/points/npts,ngroups,nhists
      common/pulses/npulse,DEF
      common/heritage/iter
      COMMON/FLAGS/FITDEAD,FIXPHASE,FITAMP
      COMMON/RUNDATA/RES,FRAMES,FNORM,IRUNNO,HISTS
      COMMON/MISSCHANNELS/MM
      common/savetime/ngo,i2pwr
      COMMON/ANS/ANS,firstgo
	common/MaxPage/n,f
     	type (qwinfo)  qw  

      P=NPTS*NGROUPS
      CALL INPUT(NHISTS,NGROUPS,NPTS,P,DATUM,SIGMA,CORR,DATT)
C	write(99,*) "back from input"
      CALL START(NGROUPS,NPTS,P)
      CALL BACK(HISTS,NGROUPS,NPTS,P,DATUM,SIGMA,CORR)
      ngo=-1
      DO 999 J=1,10
      ngo=ngo+1
      write(str,'('' CYCLE NUMBER='',i2)')ngo

      call print_log_msg("notice", TRIM(str))

        write(str,345)
345     format(' iter  test:<.02?   entropy    chitarget     chisq',
     +   '      freqsum ')
      call print_log_msg("notice", TRIM(str))

1     CALL MAXENT(NGROUPS,NPTS,P,DATUM,SIGMA,def,BASE,10,.FALSE.)
      IF (fitdead) THEN
          CALL DEADFIT(NGROUPS,NPTS,P,DATUM,SIGMA
     +                  ,CORR,DATT)
      ELSE  
          CALL MODBAK(HISTS,NGROUPS,NPTS,P,DATUM,SIGMA)
      ENDIF
      IF(fixphase) THEN
          CALL MODAMP(HISTS,NGROUPS,NPTS,P,DATUM,SIGMA)
      ELSE
          CALL MODAB(HISTS,NGROUPS,NPTS,P,DATUM,SIGMA)
      ENDIF
c      if( (iter.eq.2).and.(ngo.ge.3) ) then
c        write(99,*) ' Getting there in one iteration, so stopping...'
c        goto 997
c      endif
!      if( (iter.eq.3).and.(ngo.ge.5) ) then
!        write(99,*) ' Pretty close, but chi**2 a bit tight: stop? [Y]'
!        read(5,110) an
!        IF(AN.NE.'N'.AND.AN.NE.'n') GOTO 997
!      endif        
999   CONTINUE
997   CALL OUTPUT(P)
      CALL OUTSPEC(NGROUPS,NPTS,N,P,DATUM,F,SIGMA,datt)
!      write(99,*) ' Again? [y]'
!      read(5,110) AN
110    FORMAT(A1)
!      IF(AN.NE.'N'.AND.AN.NE.'n') GOTO 111    
!      STOP      
      return
      END
