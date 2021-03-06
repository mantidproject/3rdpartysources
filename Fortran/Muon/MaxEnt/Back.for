C
C ************************************************************
C
      SUBROUTINE BACK(HISTS,NGROUPS,NPTS,P,DATUM,SIGMA,CORR)
      INTEGER P
      REAL ZR(8192),ZI(8192),AMP(96),PH(96)
      REAL CORR(NPTS,NGROUPS),HISTS(NGROUPS)
      COMMON/MISSCHANNELS/MM
      COMMON/FLAGS/FITDEAD,FIXPHASE,FITAMP
      COMMON/RUNDATA/RES,IRUNNO,IFRAM
      COMMON/DETECT/A,B,E,c,d
      COMMON/FASE/PHASE(96),phshift
      common/sense/phi(96),TAUD(96),phases(96)
      REAL A(96),B(96),c(96),d(96),E(8192),taud
      REAL DATUM(NPTS,NGROUPS),SIGMA(NPTS,NGROUPS)
      CHARACTER*1 FITAMP
	CHARACTER*2046 str      
      DATA AMP /1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
     &          1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,48*1.0/
      DO 2 J=1,NGROUPS

      IF(HISTS(J).EQ.0) THEN
       GOTO 5
      ENDIF
      Ax=0.
      Bx=0.
      DO 3 I=1,npts
      Ax=Ax+E(I)*DATUM(I,J)/(SIGMA(I,J)*SIGMA(I,J))
3     Bx=Bx+E(I)*E(I)/(SIGMA(I,J)*SIGMA(I,J))
      SCALE=Ax/Bx
5     d(j)=scale
      DO 4 I=1,npts
      IF(sigma(i,j).gt.1e6) goto 4
      DATUM(I,J)=DATUM(I,J)-SCALE*E(I)
4     continue
2     CONTINUE
c      write(99,*) ' DATUM AND EXP...'
c      write(99,*) (DATUM(100,l),l=1,4)
c      write(99,*) (E(100)*D(L),l=1,4)
C       FIRST APPROXIMATION TO AMPLITUDE : PROP TO EXP   
      GGRO=0.
      SUM=0.
      DO 10 J=1,NGROUPS
      IF(HISTS(J).EQ.0)GOTO 10
      GGRO=GGRO+1.
      SUM=SUM+D(J)
10    CONTINUE
      SUM=SUM/GGRO
      DO 11 J=1,NGROUPS
      IF(HISTS(J).EQ.0)THEN
        AMP(J)=0.0
      ELSE
        AMP(J)=D(J)/SUM
      ENDIF
11    CONTINUE
      
      DO 1 I=1,NGROUPS
      PH(I) = phases(i)
      IF(HISTS(I).EQ.0)GOTO 1
      PHASE(I)=PH(I)/57.296
      A(I)=AMP(I)*COS(PHASE(I))
      B(I)=AMP(I)*SIN(PHASE(I))
1     CONTINUE

      call print_log_msg("notice", "PHASES AS LOADED...")
	write(str,*) (PH(I),I=1,NGROUPS)
      call print_log_msg("notice", TRIM(str))

      call print_log_msg("notice", "EXPONENTIALS AS FITTED...")
	write(str,*) (D(J),J=1,NGROUPS)
      call print_log_msg("notice", TRIM(str))

      RETURN
      END
