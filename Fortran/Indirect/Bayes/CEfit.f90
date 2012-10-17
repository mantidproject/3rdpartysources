C      PROGRAM CHUDLEY_ELLIOT_BAYES_POLY
C
      SUBROUTINE cefit(nd,X_in,Y_in,E_in,sfile,l_fn,
     1 kill,res,no,XOUT,YOUT)
      INCLUDE 'jump_inc.f90'
      integer nd, l_fn
cf2py intent(in) :: nd, l_fn                    !#points & length filename
      real X_in(1000), Y_in(1000), E_in(1000)
cf2py intent(in) :: X_in, Y_in, E_in            !x,y,e
      character*140 sfile
cf2py intent(in) :: sfile                       !sample filename
      integer kill, no
cf2py intent(out) :: kill, no
      real res(6), YOUT(mpt), XOUT(mpt)
cf2py intent(out) :: res, YOUT, XOUT           !TOTALS out
      REAL             XXX(mpt,13),SIG2(mpt)
      REAL             XM(13),YLOGPM(13),AM(13,13),SIGAM(13,13)
      REAL             BIGHES(13,13),BIGAMX(13,13)
      REAL             GRAD(13),HESS(13,13),COVAR(13,13)
      REAL             AMCHD(13)
      INTEGER          INDX(15)
      DATA             MMAX,NOUT /0,400/
C
      kill=0
      l_lpt=l_fn+4
      lptfile(1:l_lpt)=sfile(1:l_fn)//'.lpt'
      ndat=nd
      if(ndat.gt.mpt)then
       kill=1
       return
      endif
      do n=1,ndat
       X(n)=X_in(n)
       Y(n)=Y_in(n)
       E(n)=E_in(n)
      end do
      call open_f(6,lptfile)
      CALL DATINT(SIG2)
      WRITE(6,1001)NDAT
1001  format(' cefit> No. of data read in = ',i5)
      close(unit=6)
      CALL XXXINT(X,XXX,NDAT,1,XOUT,NOUT,XOFSET,IOUT1,NOUT1)
      CALL CHDLEY(SIG2,ACE,XKCE,Q,RR,XOUT,YOUT,NOUT)
      CALL HESINT(SIG2,NDAT,1,BIGHES,XXX,BIGAMX,RR,XMU0)
      DETLGA=0.0
      MORDER=0
        AMCHD(MORDER+1)=ACE
        AMCHD(MORDER+2)=XKCE
        CALL VCOPY(AMCHD,AM(1,MORDER+1),MORDER+2)
        CHIOLD=CCHI(AM(1,MORDER+1),XXX)
        DO 20 I=1,20
          CALL FILMTX(BIGHES,MMAX,HESS,MORDER+2,2)
          CALL HESGRD(HESS,SIG2,FIT,NDAT,MORDER+2,GRAD,XXX,
     *                AM(1,MORDER+1))
          CALL INVERT(HESS,COVAR,MORDER+2,INDX,DETLGB)
          CALL NEWEST(COVAR,GRAD,MORDER+2,AM(1,MORDER+1))
          CHISQ=CCHI(AM(1,MORDER+1),XXX)
          CHIDIF=CHIOLD-CHISQ
          IF (CHIDIF.LT.-0.01) THEN 
           OPEN(UNIT=6,FILE=lptfile,STATUS='old',FORM='formatted',
     1     access='append')
           write(6,1102)CHIOLD,CHISQ
1102       format(' cefit> *** Dodgy refinement :',2e13.4)
           close(unit=6)
           STOP
          ELSEIF(CHIDIF.LT.0.0001) THEN
            IF (I.GT.5) GOTO 1
          ENDIF
          CHIOLD=CHISQ
  20    CONTINUE
   1    CALL ERRBAR(COVAR,MORDER+2,SIGAM(1,MORDER+1),CHISQ)
        CALL FILMTX(BIGHES,MMAX,HESS,MORDER+2,2)
        CALL HESGRD(HESS,SIG2,FIT,NDAT,MORDER+2,GRAD,XXX,
     *              AM(1,MORDER+1))
        CALL ADPRBA(HESS,MORDER+2,BIGAMX,MMAX,XMU0)
        CALL INVERT(HESS,COVAR,MORDER+2,INDX,DETLGB)
        CALL PROBM(DETLGA,DETLGB,CHISQ,NDAT,XMU0,MORDER,XM,YLOGPM)
      CHISQ=CCHI(AM(1,1),XXX)
      CALL POLYFT(AM(1,1),0,XOUT,YOUT,NOUT,XOFSET)
      OPEN(UNIT=6,FILE=lptfile,STATUS='old',FORM='formatted',
     1 access='append')
      WRITE(6,120) CHISQ
 120  FORMAT(' Normalised Chi-squared = ',F12.3)
      WRITE(6,141) YLOGPM(1)
 141  FORMAT(' Log10[Prob(Chudley-Elliot|{Data})] = ',F10.1)
      WRITE(6,131) AM(1,1),SIGAM(1,1)
 131  FORMAT(' Coeff.  A =  ',1PE10.3,' +- ',E9.2)
      WRITE(6,132) AM(2,1),SIGAM(2,1)
 132  FORMAT(' Coeff.  K =  ',1PE10.3,' +- ',E9.2)
      CLOSE(UNIT=6)
      res(1)=CHISQ
      res(2)=YLOGPM(1)
      res(3)=AM(1,1)
      res(4)=SIGAM(1,1)
      res(5)=AM(2,1)
      res(6)=SIGAM(2,1)
      no=NOUT
      RETURN
      END
C     --------------------
      FUNCTION CCHI(A,XXX)
      INCLUDE 'jump_inc.f90'
      REAL  A(*),XXX(NDAT,*)
      CHI=0.0
      ACE=A(MORDER+1)
      XKCE=A(MORDER+2)
      do K=1,NDAT
        XK=XKCE*X(K)
        SUMm=ACE*(1.0-SIN(XK)/XK)
        do I=1,MORDER
          SUMm=SUMm+A(I)*XXX(K,I)
        end do
        FIT(K)=SUMm
        IF (E(K).GT.1.0E-10) THEN
          DIF=(FIT(K)-Y(K))/E(K)
          CHI=CHI+DIF*DIF
        ENDIF
      end do
      CCHI=CHI/FLOAT(NDAT)
      END
C
C***<chudley-elliot>***************************************************
C
      SUBROUTINE CHDLEY(SIG2,ACE,XKCE,Q,RR,XOUT,YOUT,NOUT)
      INCLUDE 'jump_inc.f90'
      REAL  SIG2(*),XOUT(*),YOUT(*)
      XMAX=X(1)
      IF (X(NDAT).GT.XMAX) XMAX=X(NDAT)
      IF (XMAX.LT.0.001)then
       OPEN(UNIT=6,FILE=lptfile,STATUS='old',FORM='formatted',
     1  access='append')
       write(6,1000)
1000   format(' cefit> something wrong here!')
       close(unit=6)
       STOP
      endif
      CMIN=1.0E+20
      NPT=500
      DXKCE=(5.0*3.141592654/XMAX)/FLOAT(NPT)
      XKCE=0.2*DXKCE
      do I=1,NPT
        XKCE=XKCE+DXKCE
        CALL ELLIOT(SIG2,XKCE,A,CHI)
        IF (CHI.LT.CMIN) THEN
          CMIN=CHI
          BESTXK=XKCE
          BESTA=A
        ENDIF
      end do
      CALL ELLIOT(SIG2,BESTXK-DXKCE,A,CM)
      C0=CMIN
      CALL ELLIOT(SIG2,BESTXK+DXKCE,A,CP)
      XKCE=BESTXK+(CM-CP)*DXKCE/(2.0*(CM-2.0*C0+CP))
      CALL ELLIOT(SIG2,XKCE,ACE,CMIN)
      IF (CMIN.GT.C0) then
       OPEN(UNIT=6,FILE=lptfile,STATUS='old',FORM='formatted',
     1  access='append')
        write(6,1102)C0,CMIN
1102    format(' cefit> Minimisation problem!! :',2e13.4)
       close(unit=6)
       STOP
      endif
      SUM1=0.0
      SUM2=0.0
      do I=1,NDAT
        XI=XKCE*X(I)
        FIT(I)=ACE*(1.0-SIN(XI)/XI)
        YDIF=Y(I)-FIT(I)
        SUM1=SUM1+YDIF
        SUM2=SUM2+YDIF*YDIF
      end do
      Q=SUM1/FLOAT(NDAT)
      RR=SUM2/FLOAT(NDAT)
      do I=1,NOUT
        XI=XKCE*XOUT(I)
        YOUT(I)=ACE*(1.0-SIN(XI)/XI)
      end do
      END
C     ---------------------------------------------
      SUBROUTINE ELLIOT(SIG2,XKCE,A,CHI)
      INCLUDE 'jump_inc.f90'
      REAL SIG2(*)
      C=0.0
      D=0.0
      do I=1,NDAT
        XI=XKCE*X(I)
        YI=1.0-SIN(XI)/XI
        YISIG2=SIG2(I)*YI
        C=C+YISIG2*Y(I)
        D=D+YISIG2*YI
      end do
      A=C/D
      CHI=A*(A*D-2.0*C)
      END
C     -----------------------------------
      SUBROUTINE POLYFT(A,M,X,Y,N,XOFSET)
      REAL A(*),X(*),Y(*)
      ACE=A(M+1)
      XKCE=A(M+2)
      do K=1,N
        XK=XKCE*X(K)
        SUMm=ACE*(1.0-SIN(XK)/XK)
        XX=1.0
        do I=1,M
          SUMm=SUMm+A(I)*XX
          XX=XX*(X(K)-XOFSET)
        end do
        Y(K)=SUMm
      end do
      END
C     -------------------------------------------------
      SUBROUTINE HESGRD(HESS,SIG2,F,N,M,GRAD,XXX,A)
      INCLUDE 'jump_inc.f90'
      REAL HESS(M,M),SIG2(*),F(*),GRAD(*),XXX(N,M),A(*)
      REAL RES(1000),WRK(1000,3)
      ACE=A(M-1)
      XKCE=A(M)
      do K=1,N
        RES(K)=SIG2(K)*(Y(K)-F(K))
        XK=XKCE*X(K)
        IF (XK.LT.1.0E-5) THEN
          WRK(K,1)=0.0
          WRK(K,2)=0.0
          WRK(K,3)=X(K)*X(K)
        ELSE
          XKK=XKCE*XK
          SXK=SIN(XK)
          CXK=COS(XK)
          WRK(K,1)=1.0-SXK/XK
          WRK(K,2)=(XK*CXK-SXK)/XKK
          WRK(K,3)=(XK*XK*SXK+2.0*(XK*CXK-SXK))/(XKCE*XKK)
        ENDIF
      end do
      do J=1,M-2
        SUMG=0.0
        SUMA=0.0
        SUMK=0.0
        do K=1,N
          SUMG=SUMG+RES(K)*XXX(K,J)
          SUMA=SUMA+SIG2(K)*XXX(K,J)*WRK(K,1)
          SUMK=SUMK+SIG2(K)*XXX(K,J)*WRK(K,2)*ACE
        end do
        GRAD(J)=-SUMG
        HESS(J,M-1)=SUMA
        HESS(M-1,J)=SUMA
        HESS(J,M)=-SUMK
        HESS(M,J)=-SUMK
      end do
      SUMA=0.0
      SUMK=0.0
      SUMAA=0.0
      SUMKK=0.0
      SUMAK=0.0
      do K=1,N
        SUMA=SUMA+RES(K)*WRK(K,1)
        SUMK=SUMK+RES(K)*WRK(K,2)*ACE
        SUMAA=SUMAA+SIG2(K)*WRK(K,1)*WRK(K,1)
        SUMKK=SUMKK+ACE*(ACE*SIG2(K)*WRK(K,2)*WRK(K,2)-RES(K)*WRK(K,3))
        SUMAK=SUMAK+WRK(K,2)*(ACE*SIG2(K)*WRK(K,1)+RES(K))
      end do
      GRAD(M-1)=-SUMA
      GRAD(M)=SUMK
      HESS(M-1,M-1)=SUMAA
      HESS(M,M)=SUMKK
      HESS(M-1,M)=-SUMAK
      HESS(M,M-1)=-SUMAK
      CALL VCOPY(HESS,WRK,M*M)
      CALL TRED2(WRK,M,M,WRK(1,2),WRK(1,3))
      CALL TQLI(WRK(1,2),WRK(1,3),M,M,WRK)
      EMAX=-1.0E20
      do I=1,M
        IF (WRK(I,2).GT.EMAX) EMAX=WRK(I,2)
      end do
      BEEFUP=EMAX/1.0E+6
      do I=1,M
        HESS(I,I)=HESS(I,I)+BEEFUP
      end do
      END
