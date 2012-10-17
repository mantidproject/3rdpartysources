      SUBROUTINE set_geom(JOM,Jann,WIDTH,WIDTH2,HEIGHT,THICK,LOUT)
      include 'muscat_com.f90'
      REAL RAD(3)
      LOGICAL LOUT
      nreg=1
      if(NREG.GT.mreg) then
       write(2,130)
       STOP
      endif
      DO J=1,mreg
       DO I=1,msurf
        IGEOM(I,J)=-1
       end do
       A(J)=0.
       B(J)=0.
       C(J)=0.
       D(J)=0.
       E(J)=0.
       F(J)=0.
       P(J)=0.
       Q(J)=0.
       R(J)=0.
      end do
      GO TO (10,20,30),JOM
10    if(LOUT)write(2,121) THICK      !INFINITE PLATE
121   FORMAT(/'   Infinite Flat Plate : Thickness = ',F8.3,/)
      G(2)=1.0E4
      G(1)=-G(2)
      G(4)=G(2)
      G(3)=-G(4)
      GOTO 25
20    if(LOUT)write(2,124)THICK,WIDTH,HEIGHT      !FINITE PLATE
124   FORMAT(/'   Finite Flat Plate: Thickness = ',F8.3
     1 ' ; Width = ',F8.3,' ; Height = ',F8.3,/)
      G(2)=0.5*HEIGHT
      G(1)=-G(2)
      G(4)=0.5*WIDTH
      G(3)=-G(4)
25    F(1)=1.
      F(2)=1.
      D(3)=1.
      D(4)=1.
      B(5)=1.
      B(6)=1.
      NSURF=4+2*NREG
      G(5)=0.
      G(6)=-THICK
      IGEOM(1,2)=1
      IGEOM(1,4)=1
      IGEOM(1,5)=1
      GOTO 35
30    if(VKz.NE.0.0) then             !CYLINDER changed from VKY??
       write(2,129)
129   FORMAT(' minus> ERROR *** PROGRAM STOPPED BECAUSE VKY ',
     1'IS NOT = 0.0'
     2  /'   WITH CYLINDRICAL SAMPLES THE PROBLEM MUST'
     3  ' BE REARRANGED SO THAT VKY = 0.0')
       STOP
      endif
      if(LOUT)then
       if(Jann.eq.0)write(2,126)
       if(Jann.eq.1)write(2,127)
126   FORMAT(/'  Cylindrical Sample ')
127   FORMAT(/'  Annular Sample ')
       write(2,128)(N,WIDTH,WIDTH2,HEIGHT,N=1,NREG)
128   FORMAT('   Region',I3,5X,'Diameter1 = ',F8.3,
     1 ' ; Diameter2 = ',F8.3,' ; Height =',F8.3,/)
      endif

C  IF NREG>1 ,REG1 IS INNERMOST CYLINDER(SAMPLE) & REG2 NEXT RING(CONT)
C   REG3 WOULD BE OUTSIDE REG2 & SO ON
      NSURF=NREG+2
      if(Jann.eq.1)NSURF=NSURF+1
      F(1)=1.
      F(2)=1.
      G(2)=0.5*HEIGHT
      G(1)=-G(2)
      IGEOM(1,2)=1
      RAD(1)=WIDTH
      RAD(2)=WIDTH2
      DO I=1,NREG
       NR1=I+2
       A(NR1)=1.
       C(NR1)=1.
       G(NR1)=-RAD(I)*RAD(I)
       G(NSURF+I)=RAD(I)
       if(Jann.eq.1)then
        NR2=I+3
        A(NR2)=1.
        C(NR2)=1.
        G(NR2)=-RAD(2)*RAD(2)
        IGEOM(1,NSURF)=1
       endif
      end do
      if(NREG.GT.1)then
       IGEOM(1,4)=0
       IGEOM(2,2)=1
       IGEOM(2,3)=1
      endif
35    continue
c      if(LOUT)then
c       write(2,132)
c132   FORMAT(/' SURFACE   x2   x    y2   y    z2   z     G',
c     1 '         GEOMETRY')
c       DO I=1,NSURF
c        write(2,131)I,A(I),B(I),C(I),D(I),E(I),F(I),G(I),
c     1 (IGEOM(N,I),N=1,NREG)
c131   FORMAT(I5,3X,6F5.1,F10.5,1X,5I3)
c       end do
c      endif

      RETURN
130   FORMAT(' minus> ERROR ** TOO MANY REGIONS')
      END
C
      SUBROUTINE TESTIN
C TESTS WETHER A PT X,Y,Z IS IN REGION IREG
C    IF ISURF.NE.0, THAT SURFACE IS NOT TESTED.
      include 'muscat_com.f90'
      IN=.TRUE.
      DO 1 I=1,NSURF
      IF(IGEOM(IREG,I).EQ.0 .OR. I.EQ.ISURF)GO TO 1
      PP=A(I)*X*X+B(I)*X+C(I)*Y*Y+D(I)*Y+E(I)*Z*Z
     + +F(I)*Z+G(I)+P(I)*X*Y+Q(I)*Y*Z+R(I)*X*Z
C  IF ON SURFACE, USE FIRST DERIVATIVE
      IF(PP.EQ.0.)PP=2.*((A(I)*X+B(I))*VX+(C(I)*Y+D(I))*VY+
     +(E(I)*Z+F(I))*VZ)+P(I)*(X*VY+Y*VX)+Q(I)*(Y*VZ+Z*VY)+
     +R(I)*(Z*VX+X*VZ)
C  IF ON SURFACE AND PARALLEL TO IT, TRY SECOND DERIVATIVE
      IF(PP.EQ.0.)PP=A(I)*VX*VX+C(I)*VY*VY+E(I)*VZ*VZ+
     + P(I)*VX*VY+Q(I)*VY*VZ+R(I)*VZ*VX
      IF(PP*IGEOM(IREG,I).LT.0)GO TO 2
    1 CONTINUE
      RETURN
   2  IN=.FALSE.
      RETURN
      END
C
      FUNCTION DIST(I)
C  FIND POSITIVE DISTANCE TO SURFACE I
C    IF ISURF.NE.0, ASSUME CURRENTLY ON SURFACE ISURF
      include 'muscat_com.f90'
C CALCULATES THE DIST TO THE I'TH SURFACE FROM (X,Y,Z) IN DIR'N
C (VX,VY,VZ)
      AA=A(I)*VX*VX+C(I)*VY*VY+E(I)*VZ*VZ+P(I)*VX*VY+Q(I)*VY*VZ
     ++R(I)*VX*VZ
      BB=2*(A(I)*VX*X+C(I)*VY*Y+E(I)*VZ*Z)
     ++B(I)*VX+D(I)*VY+F(I)*VZ+P(I)*(X*VY+Y*VX)
     ++Q(I)*(Y*VZ+Z*VY)+R(I)*(X*VZ+Z*VX)
      CC=0.
      IF(I.NE.ISURF)CC=A(I)*X*X+B(I)*X+C(I)*Y*Y+D(I)*Y+E(I)*Z*Z+
     + F(I)*Z+G(I)+P(I)*X*Y+Q(I)*Y*Z+R(I)*X*Z
C WE NOW HAVE TO SOLVE THE QUADRATIC AA*T*T+BB*T+CC=0
      IF(AA.NE.0.)GO TO 1
      IF(BB.NE.0.)GO TO 2
C  TRAJECTORY PARALLEL TO PLANE, NO INTERSECTION
      DIST=-1.0E4
      RETURN
C  SURFACE HAS NO CURVATURE IN DIRECTION OF TRAJECTORY
  2   DIST=-CC/BB
      RETURN
C  FULL QUADRATIC TREATMENT REQUIRED
1     BB=-0.5*BB/AA
      DD=BB*BB-CC/AA
      IF(DD.GE.0.)GO TO 3
C  TRAJECTORY DOES NOT INTERSECT THIS QUADRATIC SURFACE
      DIST=-2.0E4
      RETURN
3     IF (CC.EQ.0.) GO TO 6
      IF (BB.EQ.0.) GO TO 4
C  TEST SIZE OF RADICAL, COMPARED TO BB
      EE=1.-CC/AA/(BB*BB)
      IF (ABS(EE).GT.1.E-4) GO TO 4
C  BETTER PRECISION BY EXPANDING SQUARE ROOT
      DD=ABS(BB)*(1.+0.5*EE)
      GO TO 5
C  SOLVE BY QUADRATIC FORMULA
4     DD=SQRT(DD)
5      DIST=BB-DD
C  IF >0, THIS IS SMALLER POSITIVE DISTANCE; ELSE, TRY LARGER ONE
      IF (DIST.LE.0) DIST=BB+DD
      RETURN
C  PRESENTLY ON THIS SURFACE, BUT MAY HIT IT AGAIN
6     DIST=BB+ABS(BB)
      RETURN
C WE HAVE NOW RETURNED THE SMALLEST +VE DIST - IF THERE IS ONE
      END
C
      SUBROUTINE DTOEX
C CALCULATES THE DISTANCE (EXDIST) TO THE EXIT OF A REGION (IREG)
C ALONG THE DIRECTION VX,VY,VZ FROM X,Y,Z.  IF IREG = 0, FINDS
C DISTANCE TO NEAREST SURFACE.  THE EXIT SURFACE IS
C RETURNED IN ISURF.  IF ISURF .NE. 0 AT ENTRY, THEN THE POINT
C IS ASSUMED TO BE ON ISURF AND THAT DISTANCE IS 0.
      include 'muscat_com.f90'
      EXDIST=100000.0
      II=0
      DO 1 I=1,NSURF
      IF (IREG.EQ.0) GO TO 2
      IF(IGEOM(IREG,I).EQ.0) GO TO 1
2     DD=DIST(I)
      IF(DD .LE. 0.0) GO TO 1
      IF(DD.GT.EXDIST) GO TO 1
      II=I
      EXDIST=DD
  1   CONTINUE
      ISURF=II
      IF(EXDIST.EQ.100000.0) EXDIST=-1.0
      RETURN
      END
C
      SUBROUTINE WHICHR
C FINDS WHICH REGION X,Y,Z IS IN & RETURNS REGION IN IREG
C IF MORE THAN 1 REGION - IREG SET = -1
C IF NO REGION RETURNS IREG=0
      include 'muscat_com.f90'
      LOGICAL FOUND
      FOUND=.FALSE.
      IT=0
C  SET ISURF=0 SO POINT WON'T BE ASSUMED TO BE ON SURFACE
      ISURF=0
      DO 1 IREG=1,NREG
      CALL TESTIN
      IF(.NOT.IN)GOTO 1
      IF(IN.AND.FOUND) GOTO 2
      IT=IREG
      FOUND=.TRUE.
   1  CONTINUE
      IREG=IT
      RETURN
  2   IREG=-1
      RETURN
      END
C
      SUBROUTINE NEXTRG
C  FINDS WHICH REGION WILL BE ENTERED ACROSS SURFACE ISURF
C    IN DIRECTION (VX,VY,VZ) AT POINT (X,Y,Z)
      include 'muscat_com.f90'
C
C  IF DON'T HAVE ISURF, DEFAULT TO CALL WHICHR
      IF(ISURF.GT.0 .AND. ISURF.LE.NSURF) GO TO 2
      CALL WHICHR
      RETURN
C
C  FIRST TEST REGIONS ON OPPOSITE SIDE OF ISURF FROM IREG
2     IG=IGEOM(IREG,ISURF)
      DO 5 IREG=1,NREG
      IF (IG*IGEOM(IREG,ISURF).GE.0) GO TO 5
      CALL TESTIN
      IF (IN) RETURN
5     CONTINUE
C
C  DIDN'T FIND IT, SO TEST REGIONS NOT BOUNDED BY ISURF
      DO 10 IREG=1,NREG
      IF(IGEOM(IREG,ISURF).NE.0)GO TO 10
      CALL TESTIN
      IF (IN) RETURN
10    CONTINUE
C
C  STILL DIDN'T FIND IT, GIVE UP
      IREG=0
      RETURN
      END

C   HARWELL routines
      FUNCTION FA01AS(I)
      COMMON/FA01ES/G
      DOUBLE PRECISION G
      G=DMOD(G* 9228907.D0,4294967296.D0)
      IF(I.GE.0)FA01AS=G/4294967296.D0
      IF(I.LT.0)FA01AS=2.D0*G/4294967296.D0-1.D0
      RETURN
      END
      SUBROUTINE FA01BS(MAXx,NRAND)
      NRAND=INT(FA01AS(1)*FLOAT(MAXx))+1
      RETURN
      END
      SUBROUTINE FA01CS(IL,IR)
      COMMON/FA01ES/G
      DOUBLE PRECISION G
      IL=G/65536.D0
      IR=G-65536.D0*FLOAT(IL)
      RETURN
      END
      SUBROUTINE FA01DS(IL,IR)
      COMMON/FA01ES/G
      DOUBLE PRECISION G
      G=65536.D0*FLOAT(IL)+FLOAT(IR)
      RETURN
      END
      BLOCK DATA
      COMMON/FA01ES/G
      DOUBLE PRECISION G
      DATA G/1431655765.D0/
      END
