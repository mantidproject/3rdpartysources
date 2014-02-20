C
C ****************************************************************
C
      SUBROUTINE MODAMP(HISTS,NGROUPS,NPTS,P,DATUM,SIGMA)
      INTEGER P
      REAL DATUM(NPTS,NGROUPS),F(68000),SIGMA(NPTS,NGROUPS),PHASE(64)   	
      REAL CS(64),SN(64),AMP(64),phi(64)
      REAL ZR(8192),ZI(8192),HISTS(NGROUPS)
      REAL A(64),B(64),c(64),d(64),E(8192)
	character*2046 str
      COMMON/MISSCHANNELS/MM
      COMMON/DETECT/A,B,E,c,d
      COMMON/FASE/PHASE,phshift
      common/amps/amp
      common/MaxPage/n,f
      CALL ZFT(P,N,F,ZR,ZI)
      V=0.
      DO 1 I=1,NGROUPS
      IF(HISTS(I).EQ.0)GOTO 1
      S=0.0
      T=0.0
      CS(I)=COS(PHASE(I))
      SN(I)=SIN(PHASE(I))
      DO 2 J=1,npts
      HIJ=CS(I)*ZR(J)+SN(I)*ZI(J)
      S=S+DATUM(J,I)*HIJ/(SIGMA(J,I)*SIGMA(J,I))
      T=T+(HIJ*HIJ)/(SIGMA(J,I)*SIGMA(J,I))
2     CONTINUE
      AMP(I)=S/T
      V=V+SQRT(AMP(I)*AMP(I))
1     CONTINUE
      V=V/FLOAT(NGROUPS-MM)      
      DO 3 I=1,NGROUPS
      AMP(I)=AMP(I)/V
      A(I)=AMP(I)*CS(I)
      B(I)=AMP(I)*SN(I)
      phi(i)=phase(i)*57.296
3     CONTINUE

      write(str,4)
4     FORMAT(1X,'AMPLITUDES ')
      call print_log_msg("notice", TRIM(str))

      write(str,*)(AMP(I),I=1,NGROUPS)
      call print_log_msg("notice", TRIM(str))

      write(str,*) 'Fixed Phases:'
      call print_log_msg("notice", TRIM(str))

      write(str,*) (phi(I),I=1,NGROUPS)
      call print_log_msg("notice", TRIM(str))

      RETURN
      END
