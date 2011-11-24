      SUBROUTINE ACYL(MS1,ASS,ASSC,ACSC,ACC,AREAS,AREAC)
      parameter (mlam=10,mann=3,mbank=20)
      real assc(mann)
      COMMON NAN,RAD(mann+1),THETA,PI,A,B,AMUS(mann),AMUT(mann),
     1 DEN(mann),amuse(mann),amute(mann)          !add elastic values
      COMMON/PROFBL/PROFIL(50),NPROF,PRSTEP
      PRSTEP=(A-B)/(NPROF-1)
C
C     DETERMINE MAXIMUM OF THE PROFILE FUNCTION
C
      PSUM=0.
      do I=1,NPROF
       PSUM=AMAX1(PSUM,PROFIL(I))
       MS=MS1
      end do
      AREAS=0
      ASS=0
      ACC=0
      ACSC=0
      NR1=NAN-1
      IF (NR1.LT.1) THEN
        RAD1=RAD(1)
        RAD2=RAD(2)
C
C  NO STEPS ARE CHOSEN SO THAT STEP WIDTH IS THE SAME FOR ALL ANNULI
C
      MS=MS1*(RAD2-RAD1)/(RAD(2)-RAD(1))
      if(ms.lt.1) ms=1 
       CALL SUMROM(AAAA,BBBA,AREAA,PSUM,1,1,A,RAD1,RAD2,MS)
       CALL SUMROM(AAAB,BBBB,AREAB,PSUM,1,1,B,RAD1,RAD2,MS)
       AREAS1=SIGN(AREAA,A)-SIGN(AREAB,B)
       areas=areas+areas1 
       ASS=(SIGN(AAAA,A)-SIGN(AAAB,B))+ASS
       ASS=ASS/AREAS
      ELSE
       DO 15 I=1,NR1
       ASSC(i)=1.0
       IF(amus(I).EQ.0.)GO TO 15
       RAD1=RAD(I)
       RAD2=RAD(I+1)
C
C  NO STEPS ARE CHOSEN SO THAT STEP WIDTH IS THE SAME FOR ALL ANNULI
C
       MS=MS1*(RAD2-RAD1)/(RAD(2)-RAD(1))
       if(ms.lt.1) ms=1 
       CALL SUMROM(AAAA,BBBA,AREAA,PSUM,I,1,A,RAD1,RAD2,MS)
       CALL SUMROM(AAAB,BBBB,AREAB,PSUM,I,1,B,RAD1,RAD2,MS)
       AREAS1=SIGN(AREAA,A)-SIGN(AREAB,B)
       areas=areas+areas1 
       ASS=(SIGN(AAAA,A)-SIGN(AAAB,B))+ASS
       if(areas1.gt.0.) ASSC(i)=(SIGN(BBBA,A)-SIGN(BBBB,B))/areas1
  15   CONTINUE
       ASS=ASS/AREAS
       ACC=0
       ACSC=0
       RAD1=RAD(NAN)
       RAD2=RAD(NAN+1)
       MS=MS1*(RAD2-RAD1)/(RAD(2)-RAD(1))
       if(ms.lt.1) ms=1 
       CALL SUMROM(AAAA,BBBA,AREAA,PSUMA,NAN,2,A,RAD1,RAD2,MS)
       CALL SUMROM(AAAB,BBBB,AREAB,PSUMB,NAN,2,B,RAD1,RAD2,MS)
       AREAC=SIGN(AREAA,A)-SIGN(AREAB,B)
       ACSC=(SIGN(AAAA,A)-SIGN(AAAB,B))/AREAC
       ACC=(SIGN(BBBA,A)-SIGN(BBBB,B))/AREAC
      ENDIF
      RETURN
      END
 
 
      SUBROUTINE SUMROM(AAA,BBB,AREA,PSUM,NRAD,JJJ,A2,R1,R2,MS)
      parameter (mlam=10,mann=3,mbank=20)
      DIMENSION PATH(3)
      COMMON NAN,RAD(mann+1),THETA,PI,A1,B1,AMUS(mann),AMUT(mann),
     1 DEN(mann),amuse(mann),amute(mann)
      COMMON/PROFBL/PROFIL(50),NPROF,PRSTEP
      REAL LIS(mann),LSS(mann),LIST,LISN,LSST,LSSN
      NAN1=NAN+1
      A=A2
      OMGADD=0.
      IF(A.LT.0.) OMGADD=PI
      A=ABS(A)
      IF(A.GE.R2) ALIM=R2
      IF(A.LT.R2) ALIM=A
      IF(A.LE.R1) ALIM=R1
      M1=MS*((ALIM-R1)/(R2-R1)+0.0001)
      M2=MS-M1
      AAA=0.
      BBB=0.
      AREA=0.
      PSUM=0.
C
C INTEGRATE THE BEAM PROFILE FROM 0 TO ALIMIT
C
      WIDTH=A
      IF(A.GE.R2) WIDTH=R2
      ASTEP=WIDTH/40
      P1=PROBE(0.,PROFIL,NPROF,PRSTEP,A1,B1)/2
      do I=1,40
       X=I*ASTEP
       X=SIGN(X,A2)
       P2=PROBE(X,PROFIL,NPROF,PRSTEP,A1,B1)/2
       PSUM=PSUM+P1+P2
       P1=P2
      end do
      PSUM=PSUM*ASTEP
      OAD=PI-THETA
      IF(M1.EQ.0) GO TO 41
      N=M1
      RSTEP=(ALIM-R1)/M1
      RADD=-0.5*RSTEP+R1
   31 DO 30 M=1,N
      R=M*RSTEP+RADD
      NOMEG=PI*R/RSTEP
      OMEGST=PI/NOMEG
      OMEGAD=-0.5*OMEGST+OMGADD
      AREAY=R*RSTEP*OMEGST*AMUS(NRAD)
      ICOUNT=0
      SUM1=0.
      SUM2=0.
      ARSUM=0.
      I=1
   35 IF(I.GT.NOMEG) GO TO 36
      OMEGA=I*OMEGST+OMEGAD
      D=R*SIN(OMEGA)
      IF(ABS(D).GT.A) GO TO 101
C
C DETERMINE A PROFILE VALUE FOR THIS 'D' VALUE
C
      PROB=PROBE(D,PROFIL,NPROF,PRSTEP,A1,B1)
C
C CALCULATE DISTANCE INCIDENT NEUTRON PASSES THROUGH EACH ANNULUS
C
      LIST=DIST(R,RAD(1),OMEGA)
      do J=2,NAN1
       II=J-1
       LISN=DIST(R,RAD(J),OMEGA)
       LIS(II)=LISN-LIST
       LIST=LISN
      end do
C
C CALCULATE DISTANCE SCATTERED NEUTRON PASSES THROUGH EACH ANNULUS
C
      O=OMEGA+OAD
      LSST=DIST(R,RAD(1),O)
      do J=2,NAN1
       II=J-1
       LSSN=DIST(R,RAD(J),O)
       LSS(II)=LSSN-LSST
       LSST=LSSN
      end do
      ANGLE=OMEGA*180/PI
C
C CALCULATE ABSORBTION FOR PATH THROUGH ALL ANNULI,AND THROUGH INNER ANNULI
C
      NR1= NAN-1
      PATH(1)=0
      IF(NR1.LT.1) then
       PATH(1)=PATH(1)+AMUT(1)*LIS(1)+AMUTe(1)*LSS(1)
      else
       do II=1,NR1
c	split into input (I) and scattered (S) paths
        PATH(1)=PATH(1)+AMUT(II)*LIS(II)+AMUTe(II)*LSS(II)
       end do
      endif
      PATH(3)=AMUT(NAN)*(LIS(NAN)+LSS(NAN))
      PATH(2)=PATH(1)+PATH(3)
      SUM1=SUM1+EXP(-PATH(JJJ))*PROB
      SUM2=SUM2+EXP(-PATH(JJJ+1))*PROB
      ARSUM=ARSUM+PROB
      I=I+1
      GO TO 35
  101 I=NOMEG-I+2
      GO TO 35
   36 AAA=SUM1*AREAY+AAA
      BBB=SUM2*AREAY+BBB
      AREA=AREA+ARSUM*AREAY
   30 CONTINUE
   41 IF(M2.EQ.0) RETURN
      N=M2
      RSTEP=(R2-ALIM)/M2
      RADD=-0.5*RSTEP+ALIM
      M2=0
      GO TO 31
      END

      FUNCTION DIST(R1,RAD,OMEGA)
      R=R1
      DIST=0.
      B=R*SIN(OMEGA)
      IF(ABS(B).GT.RAD) RETURN
      T=R*COS(OMEGA)
      C=RAD*RAD-B*B
      D=SQRT(C)
      IF(R.GT.RAD) GO TO 1
      DIST=T+D
      RETURN
    1 DIST=D*(1+SIGN(1.,T))
      RETURN
      END

      FUNCTION PROBE(X,PROFIL,NPROF,PRSTEP,A1,B1)
      DIMENSION PROFIL(NPROF)
      IF(X.GT.A1.OR.X.LT.B1) GO TO 10
      DIFF=A1-X
      APOS=DIFF/PRSTEP
      NPOS=INT(APOS)
      DIFF=APOS-FLOAT(NPOS)
      NPOS1=NPOS+1
      NPOS2=NPOS+2
      PROBE=PROFIL(NPOS1)+(PROFIL(NPOS2)-PROFIL(NPOS1))*DIFF
      RETURN
   10 PROBE=0.
      RETURN
      END

      FUNCTION EQIVAR(RAD,A,B)
      COMMON/PROFBL/PROFIL(50),NPROF,PRSTEP
      PRSTEP=(A-B)/(NPROF-1)
      NSTEP=2000
      STEP=2.*RAD/NSTEP
      STEPH=0.5*STEP
      SUM=0.
      DO 10 I=1,NSTEP
      X=RAD-(I-1)*STEP+STEPH
      PROB=PROBE(X,PROFIL,NPROF,PRSTEP,A,B)
      Y=RAD*RAD-X*X
      IF(Y.LE.0.) GO TO 10
      Y=2.*SQRT(Y)
      SUM=SUM+Y*PROB
   10 CONTINUE
      EQIVAR=SUM*STEP
      RETURN
      END

      FUNCTION AINTER(ALIN,COR,AL,NPTS,NACORD)
      DIMENSION ALIN(NACORD),COR(nacord)
C  Interpolate COR at the point X=AL, using linear interpolation.
      X1=ALIN(1)
      Y1=COR(1)
      DO 10 I=2,NPTS
      X1S=X1
      Y1S=Y1
      X2=ALIN(I)
      Y2=COR(I)
      IF(AL.LT.X2) GO TO 20
      X1=X2
      Y1=Y2
   10 CONTINUE
   20 CONTINUE
      GRAD=(Y2-Y1S)/(X2-X1S)
      AINTER=Y1S+GRAD*(AL-X1S)
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
