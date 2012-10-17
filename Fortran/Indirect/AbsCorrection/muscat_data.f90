      SUBROUTINE open_f(nit,fname)
      character*140 fname
      logical found   
      INQUIRE(FILE=fname,EXIST=found)
      if(found)then
       OPEN(UNIT=nit,FILE=fname,STATUS='OLD',FORM='FORMATTED')
       CLOSE(UNIT=nit,STATUS='DELETE')
      endif
      OPEN(UNIT=nit,FILE=fname,STATUS='NEW',FORM='FORMATTED')
      RETURN
      END

      SUBROUTINE Q_dir(Ukx,Uky,Ukz,QSS)
      include 'muscat_com.f90'
      parameter (pi=3.14159265, AH=0.4818025)
      M7 = 0
1008  M7 = M7+1
      if (M7.GT.1000) then
       RETURN
      endif
      Qv = Qran*FA01AS(1) + QMIN
      CALL FA01BS(NW,IRW)             !NW is full +/- range with elastic at Nel
      W = (IRW-Nel)*DW
      Vkm2 = Vkm*Vkm
      TKM2 = (Vkm2 + AH*W)
      IF(TKM2.lt.0.0)goto 1008
      TKM = SQRT(TKM2)
      COST = (TKM2 + Vkm2 - Qv*Qv)/(2.0*Vkm*TKM)
      if(COST.gt.1.0.OR.COST.lt.-1.0)goto 1008   !MOD(Cos) must be <1
      IIW = IRW                       !use full range
      CALL SINT
c      IF(W.LT.0.0) SQW=SQW*EXP(-HOVERT*W)     !why this?
      SQW=SQW*EXP(-HOVERT*W)
      QSS = Qv*SQW
      B9 = B9*SQW*Qv
      phi = FA01AS(1)*2.0*PI
      QT = TKM*SQRT(1-COST*COST)
      QK = TKM*COST - Vkm
      B2 = (QK + Vkm)/TKM
      B3 = QT/TKM
      if(UKZ.LT.1.0) then
       A2 = SQRT(1.0 - UKZ*UKZ)
       UQTZ = COS(phi)*A2
       UQTX = -COS(phi)*UKZ*UKX/A2 + SIN(phi)*UKY/A2
       UQTY = -COS(phi)*UKZ*UKY/A2 - SIN(phi)*UKX/A2
       UKX = B2*UKX + B3*UQTX
       UKY = B2*UKY + B3*UQTY
       UKZ = B2*UKZ + B3*UQTZ
      else
       UKX = B3*COS(phi)
       UKY = B3*SIN(phi)
       UKZ = B2
      endif
      Vkm = TKM
      RETURN
      END

      SUBROUTINE NEWV
      include 'muscat_com.f90'
      SIGA1 = SIGA*VKINC/VKM
      SIGS1 = SIGB
      SIGT = SIGS1 + SIGA1            !new total x-sec
      Vmu = DENS*SIGT             !new trans coeff
      VMFP = 1.0/Vmu              !mean free path
      RETURN
      END

      SUBROUTINE INIT_SQW
      include 'muscat_com.f90'
      DO I=1,mq              !initialise S(Q,w) array
       DO J=1,mw1
        S(I,J)=0.
       end do
      end do
      RETURN
      END

      SUBROUTINE LOG_SIG
      include 'muscat_com.f90'
      DO J=1,NW               !check values of S(Q,w)
       DO I=1,NQ
        if(S(I,J).GT.0.0)then
         S(I,J) = ALOG(S(I,J))
        else
         S(I,J) = -20.0
        endif
       end do
      end do
      RETURN
      END

      SUBROUTINE unit_vector(Vkx,Vky,Vkz)
      ONE = SQRT(VKX*VKX + VKY*VKY + VKZ*VKZ) !unit vector 
      VKX = VKX/ONE               !direction
      VKY = VKY/ONE
      VKZ = VKZ/ONE
      RETURN
      END
  
      SUBROUTINE SINT
      include 'muscat_com.f90'
       if(Qv.LE.Qq(1))then          !use lowest value I=1
        SQW = S(1,IIW)
       else
        if(Qv.GT.Qq(NQ))then          !if > max value
         SQW = 0.0             !set=0
        else
         IQ0 = INT((Qv-Qq(1))/DQ) + 1
         IF(IQ0.GT.NQ-2) IQ0 = NQ-2        !interpolate value
         IQ1 = IQ0+1
         IQ2 = IQ0+2
         U = (Qv-(IQ0-1)*DQ - Qq(1))/DQ
         Aa = (S(IQ0,IIW) -2*S(IQ1,IIW) +S(IQ2,IIW) )/2
         Bb = (-3*S(IQ0,IIW) +4*S(IQ1,IIW) -S(IQ2,IIW) )/2
         Cc = S(IQ0,IIW)
         SQW = EXP(Aa*U*U + Bb*U + Cc)
        endif
       endif
      RETURN
      END

      SUBROUTINE RONE(NOINC,JOM)          !get starting point
      include 'muscat_com.f90'
      Vx=Vkx
      Vy=Vky
      Vz=Vkz
      IREG=NREG
      GO TO (1000,2000,3000),JOM      !check Geom type
C
C     INFINITE FLAT PLATE SAMPLES
C
1000  X=0.                  !entry flat face
      Y=0.                  !centre  width
      Z=0.                  !centre  height
      GOTO 4000
C
C     FINITE FLAT PLATE SAMPLES
C
 2000 X = 0.0                   !entry flat face
      Y = FA01AS(-1)*G(4)           !random width
      Z = FA01AS(-1)*G(2)           !random height
      GOTO 4000
C
C     CYLINDRICAL SAMPLES
C
 3000 Y = FA01AS(-1)*G(NSURF+NREG)      !random width
      Z = FA01AS(-1)*G(2)           !random height
      X = -SQRT(-G(NSURF) - Y*Y)        !entry curved face
C
4000  isurf=0
      CALL DTOEX          !find distance to exit =exdist
      DL=exdist
      B9=1.-EXP(-Vmu*DL)          !atten to exit
      VL=-(VMFP*ALOG(1-FA01AS(1)*B9))
c   VL=FA01AS(1)*DL             !random point along DL =VL
      B9=B9/SIGT              !new weighting
      CALL inc_xyz(VL)            !new xyz position
      NOINC=NOINC+1               !increment ntn number
      RETURN
      END

      SUBROUTINE inc_xyz(VL)
      include 'muscat_com.f90'
      X = X + VL*Vx
      Y = Y + VL*Vy
      Z = Z + VL*Vz
      RETURN
      END
