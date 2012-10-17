      SUBROUTINE erange(Ndat,Xin,Yin,Ein,er,Nbin,Rscl,
     1 Nout,BNORM,Xout)
      INCLUDE 'res_par.f90'
      real, dimension(m_d),intent(in) :: Xin, Yin, Ein
      integer, intent(in) :: Ndat, Nbin
      real, dimension(2),intent(in) :: er
      real, intent(in) :: Rscl
      integer, dimension(3),intent(out) :: Nout
      real, intent(out) :: BNORM
      real, dimension(m_d),intent(out) :: Xout
      real*4 XD(m_d),D(m_d),E(m_d)
C
      EMIN=er(1)
      EMAX=er(2)
      SMALL=1.0E-20
      DO I=1,NDAT
       XD(I)=Xin(I)
       D(I)=Yin(I)
       E(I)=Ein(I)
       IF (E(I).GT.SMALL) THEN
        E(I)=E(I)*E(I)
       ELSE
        E(I)=0.0
       ENDIF
      END DO
      DO I=1,NDAT
       IF (XD(I).GE.EMIN) GOTO 2
      END DO
   2  IMIN=I
      DO I=NDAT,1,-1
       IF (XD(I).LE.EMAX) GOTO 3
      END DO
   3  IMAX=I
      BNORM=1.0/FLOAT(NBIN)
      N=0
      do I=IMIN,IMAX,NBIN
       N=N+1
       XXD=0.0
       DD=0.0
       EE=0.0
       K=0
       do J=0,NBIN-1
        XXD=XXD+XD(I+J)
        IF (E(I+J).GT.SMALL) THEN
         K=K+1
         DD=DD+D(I+J)
         EE=EE+E(I+J)
        ENDIF
       end do
       XD(N)=BNORM*XXD
       IF (K.GT.0) THEN
        D(N)=BNORM*DD
        E(N)=2.0*FLOAT(K*K)/(EE*RSCL)
       ELSE
        D(N)=0.0
        E(N)=0.0
       ENDIF
       Xout(N)=XD(N)
      end do
      Nout(1)=N
      Nout(2)=IMIN
      Nout(3)=IMAX
      RETURN
      END
