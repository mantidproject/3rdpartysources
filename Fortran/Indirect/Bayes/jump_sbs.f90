C     --------------------------------------------------------
      SUBROUTINE PROBM(DETLGA,DETLGB,CHISQ,N,XMU0,M,XM,YLOGPM)
      REAL XM(*),YLOGPM(*)
      XM(M+1)=FLOAT(M-1)
      XMU=XMU0*FLOAT(M)
      XMULG=0.0
      IF (M.GT.0) XMULG=FLOAT(M)*LOG10(XMU)/2.0
      DETLGB=DETLGB+FLOAT(M+2)*LOG10(0.25)
      CHI2LG=LOG10(2.71828185)*CHISQ*FLOAT(N)/2.0
      YLOGPM(M+1)=XMULG-CHI2LG+0.5*(DETLGA-DETLGB)
      END
C
C***<read in & initialise data>*****************************************
C
      SUBROUTINE DATINT(SIG2)
      INCLUDE 'jump_inc.f90'
      REAL  SIG2(*)
      SMALL=1.0E-10
      ESCL=1.0
      CALL VRFILL(SIG2,0.0,NDAT)
      do I=1,NDAT
        E(I)=E(I)*ESCL
        IF (E(I).GT.SMALL) SIG2(I)=2.0/(E(I)*E(I))
      end do
      END
c
C***<Gradients and things>**********************************************
C
      SUBROUTINE XXXINT(X,XXX,N,M,XOUT,NOUT,XOFSET,IOUT1,NOUT1)
      REAL X(*),XXX(N,M),XOUT(*)
C
      XMIN=X(1)
      XMAX=X(1)
      do I=2,N
        IF (X(I).LT.XMIN) XMIN=X(I)
        IF (X(I).GT.XMAX) XMAX=X(I)
      end do
      XOFSET=0.5*(XMIN+XMAX)
      do K=1,N
        XXX(K,1)=1.0
        do J=2,M
          XXX(K,J)=XXX(K,J-1)*(X(K)-XOFSET)
        end do
      end do
      DX=XMAX*1.05/FLOAT(NOUT)
      XOUT(1)=DX
      do I=2,NOUT
        XOUT(I)=XOUT(I-1)+DX
      end do
      do I=1,NOUT
        IF (XOUT(I).GE.XMIN) GOTO 1
      end do
   1  IOUT1=I
      do I=1,NOUT
        IF (XOUT(IOUT1+I-1).GE.XMAX) GOTO 2
      end do
   2  NOUT1=I
      END
C
C     --------------------------------------------------
      SUBROUTINE HESINT(SIG2,N,M,HESS,XXX,AMTRX,RR,XMU0)
      REAL SIG2(*),HESS(M,M),XXX(N,M)
      REAL AMTRX(M,M)
      REAL W(1000),D(35),E(35)
      do I=1,M
        do J=I,M
          SUM1=0.0
          SUM2=0.0
          do K=1,N
            XKJ=XXX(K,I)*XXX(K,J)
            SUM1=SUM1+SIG2(K)*XKJ
            SUM2=SUM2+XKJ
          end do
          HESS(I,J)=SUM1
          HESS(J,I)=SUM1
          AMTRX(I,J)=SUM2
          AMTRX(J,I)=SUM2
        end do
      end do
      XMU0=RR*FLOAT(N)
      XMU0=0.5/XMU0
      CALL VCOPY(HESS,W,M*M)
      CALL TRED2(W,M,M,D,E)
      CALL TQLI(D,E,M,M)
      EMAX=-1.0E20
      do I=1,M
        IF (D(I).GT.EMAX) EMAX=D(I)
      end do
      BEEFUP=EMAX/1.0E+6
      do I=1,M
        HESS(I,I)=HESS(I,I)+BEEFUP
      end do
      END
C     --------------------------------------
      SUBROUTINE FILMTX(BIGA,MMAX,A,M,IFLAG)
      REAL BIGA(MMAX,*),A(M,*)
      MFILL=M
      IF (IFLAG.GT.1) MFILL=M-2
      do J=1,MFILL
        do I=1,MFILL
          A(I,J)=BIGA(I,J)
        end do
      end do
      END
C     ------------------------------------------
      SUBROUTINE ADPRBA(HESS,M2,BIGAMX,MMAX,XMU0)
      REAL HESS(M2,M2),BIGAMX(MMAX,MMAX)
      ASCL=4.0*XMU0
      do J=1,M2-2
        do I=1,M2-2
          HESS(I,J)=HESS(I,J)+ASCL*BIGAMX(I,J)
        end do
      end do
      END
C     --------------------------------------------
      SUBROUTINE INVERT(HESS,COVAR,NP,INDX,DETLOG)
      REAL    HESS(NP,*),COVAR(NP,*)
      INTEGER INDX(*)
      SMALL=1.E-20
      DETLOG=0.0
      CALL VRFILL(COVAR,0.0,NP*NP)
      do I=1,NP
        COVAR(I,I)=1.0
      end do
      CALL LUDCMP(HESS,NP,NP,INDX,D)
      do I=1,NP
        DETLOG=DETLOG+LOG10(ABS(HESS(I,I))+SMALL)
        CALL LUBKSB(HESS,NP,NP,INDX,COVAR(1,I))
      end do
      END
C     --------------------------------------
      SUBROUTINE ERRBAR(COVAR,NP,SIGA,CHISQ)
      REAL COVAR(NP,*),SIGA(*)
      SMALL=1.0E-20
      SIGSCL=SQRT(CHISQ)
      do I=1,NP
        SIGA(I)=SIGSCL*SQRT(2.0*ABS(COVAR(I,I))+SMALL)
      end do
      END
C     -----------------------------------
      SUBROUTINE NEWEST(COVAR,GRAD,NP,EX)
      REAL  COVAR(NP,*),GRAD(*),EX(*),DV(13)
      CALL MLTMXV(GRAD,COVAR,NP,DV)
      do I=1,NP
        EX(I)=EX(I)-DV(I)
      end do
      END
C     ---------------------------
      SUBROUTINE MLTMXV(P,OP,N,D)
      REAL P(*),OP(N,N),D(*)
      do K=1,N
        SUMm=0.0
        do J=1,N
          SUMm=SUMm+OP(J,K)*P(J)
        end do
        D(K)=SUMm
      end do
      END
C
C***<numerical recipes routines>***************************************
C
      SUBROUTINE TRED2(A,N,NP,D,E)
      DIMENSION A(NP,NP),D(NP),E(NP)
      IF(N.GT.1)THEN
        do I=N,2,-1  
          L=I-1
          H=0.
          SCAL=0.
          IF(L.GT.1)THEN
            do K=1,L
              SCAL=SCAL+ABS(A(I,K))
            end do
            IF(SCAL.EQ.0.)THEN
              E(I)=A(I,L)
            ELSE
              do K=1,L
                A(I,K)=A(I,K)/SCALE
                H=H+A(I,K)**2
              end do
              F=A(I,L)
              G=-SIGN(SQRT(H),F)
              E(I)=SCAL*G
              H=H-F*G
              A(I,L)=F-G
              F=0.
              do J=1,L
                G=0.
                do K=1,J
                  G=G+A(J,K)*A(I,K)
                end do
                IF(L.GT.J)THEN
                  do K=J+1,L
                    G=G+A(K,J)*A(I,K)
                  end do
                ENDIF
                E(J)=G/H
                F=F+E(J)*A(I,J)
              end do
              HH=F/(H+H)
              do J=1,L
                F=A(I,J)
                G=E(J)-HH*F
                E(J)=G
                do K=1,J
                  A(J,K)=A(J,K)-F*E(K)-G*A(I,K)
                end do
              end do
            ENDIF
          ELSE
            E(I)=A(I,L)
          ENDIF
          D(I)=H
        end do
      ENDIF
      E(1)=0.
      do I=1,N
        D(I)=A(I,I)
      end do
      END
C     ---------------------------
      SUBROUTINE TQLI(D,E,N,NP)
      DIMENSION D(NP),E(NP)
      IF (N.GT.1) THEN
        do I=2,N
          E(I-1)=E(I)
        end do
        E(N)=0.
        do L=1,N
          ITER=0
1         do M=L,N-1
            DD=ABS(D(M))+ABS(D(M+1))
            IF (ABS(E(M))+DD.EQ.DD) GO TO 2
          end do
          M=N
2         IF(M.NE.L)THEN
            IF(ITER.EQ.30)STOP 'too many iterations'
            ITER=ITER+1
            G=(D(L+1)-D(L))/(2.*E(L))
            R=SQRT(G**2+1.)
            G=D(M)-D(L)+E(L)/(G+SIGN(R,G))
            S=1.
            C=1.
            P=0.
            do I=M-1,L,-1
              F=S*E(I)
              B=C*E(I)
              IF(ABS(F).GE.ABS(G))THEN
                C=G/F
                R=SQRT(C**2+1.)
                E(I+1)=F*R
                S=1./R
                C=C*S
              ELSE
                S=F/G
                R=SQRT(S**2+1.)
                E(I+1)=G*R
                C=1./R  
                S=S*C
              ENDIF
              G=D(I+1)-P
              R=(D(I)-G)*S+2.*C*B
              P=S*R
              D(I+1)=G+P
              G=C*R-B
            end do
            D(L)=D(L)-P
            E(L)=G
            E(M)=0.
            GO TO 1
          ENDIF
        end do
      ENDIF
      END
